package gridss;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.concurrent.ExecutorService;

import com.google.common.io.Files;

import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.FileSystemContext;
import au.edu.wehi.idsv.SingleReadEvidence;
import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;
import gridss.cmdline.VcfTransformCommandLineProgram;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Log;
import picard.cmdline.CommandLineProgramProperties;

@CommandLineProgramProperties(
        usage = "Annotates breakpoint variant calls",  
        usageShort = "Annotates breakpoint variant calls"
)
public class AnnotateVariants extends VcfTransformCommandLineProgram {
	private static final Log log = Log.getInstance(AnnotateVariants.class);
	public static void writeAssemblyBreakends(File file, AssemblyEvidenceSource assemblyEvidence) throws IOException {
		log.info("Writing breakend assembly support.");
		File tmp = gridss.Defaults.OUTPUT_TO_TEMP_FILE ? FileSystemContext.getWorkingFileFor(file) : file;
		try (BufferedOutputStream writer = new BufferedOutputStream(new FileOutputStream(tmp))) {
			try (CloseableIterator<DirectedEvidence> it = assemblyEvidence.iterator()) {
				while (it.hasNext()) {
					SingleReadEvidence ass = (SingleReadEvidence) it.next();
					if (!ass.getSAMRecord().isSecondaryOrSupplementary()) {
						writer.write('>');
						writer.write(ass.getEvidenceID().getBytes(StandardCharsets.US_ASCII));
						writer.write('\n');
						if (ass.getBreakendSummary().direction == BreakendDirection.Forward) {
							writer.write(ass.getAnchorSequence(), 0, ass.getAnchorSequence().length);
							writer.write(ass.getBreakendSequence(), 0, ass.getBreakendSequence().length);
						} else {
							writer.write(ass.getBreakendSequence(), 0, ass.getBreakendSequence().length);
							writer.write(ass.getAnchorSequence(), 0, ass.getAnchorSequence().length);
						}
						writer.write('\n');
					}
				}
			}
		}
		if (tmp != file) {
			Files.move(tmp, file);
		}
		log.info("Writing breakend assembly support complete.");
	}
	@Override
	public CloseableIterator<VariantContextDirectedBreakpoint> iterator(CloseableIterator<VariantContextDirectedBreakpoint> calls, ExecutorService threadpool) {
		AllocateEvidence ae = new AllocateEvidence();
		AnnotateReferenceCoverage arc = new AnnotateReferenceCoverage();
		AnnotateInexactHomology ihom = new AnnotateInexactHomology(); 
		copyInputs(ae);
		copyInputs(arc);
		copyInputs(ihom);
		ae.INPUT_VCF = INPUT_VCF; // needed for caching 
		calls = ae.iterator(calls, threadpool);
		calls = arc.iterator(calls, threadpool);
		calls = ihom.iterator(calls, threadpool);
		return calls;
	}
	public static void main(String[] argv) {
        System.exit(new AnnotateVariants().instanceMain(argv));
    }
}
