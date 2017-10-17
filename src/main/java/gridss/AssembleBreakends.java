package gridss;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.ExecutorService;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;

import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.SAMEvidenceSource;
import gridss.cmdline.MultipleSamFileCommandLineProgram;
import htsjdk.samtools.util.IOUtil;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(
		summary = "Assembles breakend contigs using positional de Bruijn graph assembly of "
        		+ "reads supporting putative structural variants.",
        oneLineSummary = "Assembles breakend contigs.",
        programGroup = gridss.cmdline.programgroups.Assembly.class
)
public class AssembleBreakends extends MultipleSamFileCommandLineProgram {
    @Argument(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output file containing subset of input", optional=false)
    public File OUTPUT;
	public static void main(String[] argv) {
        System.exit(new AssembleBreakends().instanceMain(argv));
    }
	@Override
	public int doWork(ExecutorService threadpool) throws IOException {
		IOUtil.assertFileIsWritable(OUTPUT);
		ProcessingContext pc = getContext();
		List<SAMEvidenceSource> sources = getSamEvidenceSources();
    	AssemblyEvidenceSource assembler = new AssemblyEvidenceSource(pc, sources, OUTPUT);
    	assembler.assembleBreakends(threadpool);
    	return 0;
	}
}
