package gridss;

import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.SAMEvidenceSource;
import gridss.cmdline.MultipleSamFileCommandLineProgram;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.ExecutorService;

@CommandLineProgramProperties(
		summary = "Assembles breakend contigs using positional de Bruijn graph assembly of "
        		+ "reads supporting putative structural variants.",
        oneLineSummary = "Assembles breakend contigs.",
        programGroup = gridss.cmdline.programgroups.Assembly.class
)
public class AssembleBreakends extends MultipleSamFileCommandLineProgram {
    @Argument(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output file containing subset of input", optional=false)
    public File OUTPUT;
	@Argument(doc="Used for scaling assembly across multiple jobs. This is the zero-based index of this job.", optional=true)
    public int JOB_INDEX = 0;
	@Argument(doc="Used for scaling assembly across multiple jobs. " +
			"This is the total number of jobs to spread the assembly over. " +
			"Work will be allocated across all jobs based on an even distribution of genomic regions to process. " +
			"After all jobs have completed, output should be gathered by rerunning AssembleBreakends with JOB_NODES=1.", optional=true)
	public int JOB_NODES = 1;
	public static void main(String[] argv) {
        System.exit(new AssembleBreakends().instanceMain(argv));
    }
	@Override
	public int doWork(ExecutorService threadpool) throws IOException {
		IOUtil.assertFileIsWritable(OUTPUT);
		ProcessingContext pc = getContext();
		List<SAMEvidenceSource> sources = getSamEvidenceSources();
    	AssemblyEvidenceSource assembler = new AssemblyEvidenceSource(pc, sources, OUTPUT);
    	assembler.assembleBreakends(threadpool, JOB_INDEX, JOB_NODES);
    	return 0;
	}
	@Override
	protected String[] customCommandLineValidation() {
		String[] val = multiSamFileInputCustomCommandLineValidation();
		if (val != null) return val;
		String[] err = getWorkingDirectoryFilenameCollisions(OUTPUT, "OUTPUT");
		if (err != null) {
			return err;
		}
		if (JOB_NODES < 1) {
			return new String[] { "JOB_NODES must be at least 1."};
		}
		if (JOB_INDEX < 0) {
			return new String[] { "JOB_INDEX cannot be less than 0."};
		}
		if (JOB_INDEX >= JOB_NODES) {
			return new String[] { "JOB_INDEX is zero-based: JOB_INDEX must be less than JOB_NODES."};
		}
		return super.customCommandLineValidation();
	}
}
