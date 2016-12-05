package gridss;

import java.io.File;
import java.io.IOException;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;

import au.edu.wehi.idsv.VariantCaller;
import gridss.cmdline.FullEvidenceCommandLineProgram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(
        usage = "Identifies putative structural variants.",  
        usageShort = "Identifies putative structural variants."
)
public class IdentifyVariants extends FullEvidenceCommandLineProgram {
	@Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="VCF structural variation calls.")
    public File OUTPUT_VCF;
	private void callVariants(ExecutorService threadpool) throws IOException, InterruptedException, ExecutionException {
		File rawCalls = getContext().getFileSystemContext().getBreakpointVcf(OUTPUT_VCF);
		if (!rawCalls.exists()) {
			VariantCaller caller = new VariantCaller(getContext(), getSamEvidenceSources(), getAssemblySource());
			caller.callBreakends(rawCalls, threadpool);
		}
	}
	public static void main(String[] argv) {
        System.exit(new IdentifyVariants().instanceMain(argv));
    }
	@Override
	public int doWork(ExecutorService threadpool) throws IOException, InterruptedException, ExecutionException {
		IOUtil.assertFileIsWritable(OUTPUT_VCF);
    	callVariants(threadpool);
		return 0;
	}
}
