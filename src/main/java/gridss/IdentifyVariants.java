package gridss;

import au.edu.wehi.idsv.VariantCaller;
import gridss.cmdline.FullEvidenceCommandLineProgram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.io.IOException;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;

@CommandLineProgramProperties(
		summary = "Identifies putative structural variants.",  
        oneLineSummary = "Identifies putative structural variants.",
        programGroup = gridss.cmdline.programgroups.VariantCalling.class
)
public class IdentifyVariants extends FullEvidenceCommandLineProgram {
	private static final Log log = Log.getInstance(IdentifyVariants.class);
	@Argument(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="VCF structural variation calls.")
    public File OUTPUT_VCF;
	@Argument(doc="Used for scaling assembly across multiple jobs. This is the zero-based index of this job.", optional=true)
	public int JOB_INDEX = 0;
	@Argument(doc="Used for scaling assembly across multiple jobs. " +
			"This is the total number of jobs to spread the assembly over. " +
			"Work will be allocated across all jobs based on an even distribution of genomic regions to process. " +
			"After all jobs have completed, output should be gathered by rerunning AssembleBreakends with JOB_NODES=1.", optional=true)
	public int JOB_NODES = 1;

	public static void main(String[] argv) {
        System.exit(new IdentifyVariants().instanceMain(argv));
    }
	@Override
	public int doWork(ExecutorService threadpool) throws IOException, InterruptedException, ExecutionException {
		IOUtil.assertFileIsWritable(OUTPUT_VCF);
		VariantCaller caller = new VariantCaller(getContext(), getSamEvidenceSources(), getAssemblySource());
		caller.callBreakends(OUTPUT_VCF, threadpool, JOB_INDEX, JOB_NODES);
		log.info("Raw variant calls written to " + OUTPUT_VCF);
		return 0;
	}
}
