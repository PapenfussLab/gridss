package gridss;

import au.edu.wehi.idsv.*;
import com.google.common.io.Files;
import gridss.cmdline.VcfTransformCommandLineProgram;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.concurrent.ExecutorService;

@CommandLineProgramProperties(
        summary = "Annotates breakpoint variant calls",  
        oneLineSummary = "Annotates breakpoint variant calls",
        programGroup = gridss.cmdline.programgroups.VariantCalling.class
)
public class AnnotateVariants extends VcfTransformCommandLineProgram {
	private static final Log log = Log.getInstance(AnnotateVariants.class);
	public static void writeAssemblyBreakends(File file, AssemblyEvidenceSource assemblyEvidence) throws IOException {
		log.info("Writing breakend assembly support.");
		File tmp = gridss.Defaults.OUTPUT_TO_TEMP_FILE ? FileSystemContext.getWorkingFileFor(file) : file;
		try (BufferedOutputStream writer = new BufferedOutputStream(new FileOutputStream(tmp))) {
			try (CloseableIterator<DirectedEvidence> it = assemblyEvidence.iterator(SAMEvidenceSource.EvidenceSortOrder.EvidenceStartPosition)) {
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
	public CloseableIterator<VariantContextDirectedEvidence> iterator(CloseableIterator<VariantContextDirectedEvidence> calls, ExecutorService threadpool) {
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
