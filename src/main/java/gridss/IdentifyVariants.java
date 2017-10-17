package gridss;

import java.io.File;
import java.io.IOException;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;

import au.edu.wehi.idsv.VariantCaller;
import gridss.cmdline.FullEvidenceCommandLineProgram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(
		summary = "Identifies putative structural variants.",  
        oneLineSummary = "Identifies putative structural variants.",
        programGroup = gridss.cmdline.programgroups.VariantCalling.class
)
public class IdentifyVariants extends FullEvidenceCommandLineProgram {
	private static final Log log = Log.getInstance(IdentifyVariants.class);
	@Argument(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="VCF structural variation calls.")
    public File OUTPUT_VCF;
	public static void main(String[] argv) {
        System.exit(new IdentifyVariants().instanceMain(argv));
    }
	@Override
	public int doWork(ExecutorService threadpool) throws IOException, InterruptedException, ExecutionException {
		IOUtil.assertFileIsWritable(OUTPUT_VCF);
		VariantCaller caller = new VariantCaller(getContext(), getSamEvidenceSources(), getAssemblySource());
		caller.callBreakends(OUTPUT_VCF, threadpool);
		log.info("Raw variant calls written to " + OUTPUT_VCF);
		return 0;
	}
}
