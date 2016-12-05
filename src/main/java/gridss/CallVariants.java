package gridss;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import com.google.common.base.Function;
import com.google.common.collect.Lists;

import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.SAMEvidenceSource;
import au.edu.wehi.idsv.util.FileHelper;
import gridss.cmdline.FullEvidenceCommandLineProgram;
import gridss.cmdline.MultipleSamFileCommandLineProgram;
import htsjdk.samtools.SamPairUtil.PairOrientation;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import picard.analysis.InsertSizeMetrics;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

/**
 * Extracts structural variation evidence and assembles breakends
 * @author Daniel Cameron
 *
 */
@CommandLineProgramProperties(
        usage = "Calls structural variations from one or more SAM/BAM input files.",  
        usageShort = "Calls structural variations from NGS sequencing data"
)
public class CallVariants extends FullEvidenceCommandLineProgram {
	private static final Log log = Log.getInstance(CallVariants.class);
	@Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="VCF structural variation calls.")
    public File OUTPUT;
	public CallVariants() {
		super(false);
	}
	private void extractEvidence(ExecutorService threadpool, List<SAMEvidenceSource> samEvidence) throws InterruptedException, ExecutionException {
		log.info("Extracting evidence.");
		for (Future<Void> future : threadpool.invokeAll(Lists.transform(samEvidence, new Function<SAMEvidenceSource, Callable<Void>>() {
				@Override
				public Callable<Void> apply(final SAMEvidenceSource input) {
					return new Callable<Void>() {
						@Override
						public Void call() throws Exception {
							try {
								input.ensureMetrics();
								InsertSizeMetrics ism = input.getMetrics().getInsertSizeMetrics();
					    		if (ism != null && ism.PAIR_ORIENTATION != PairOrientation.FR) {
					    			String msg = "GRIDSS currently supports only FR read pair orientation. If usage with other read pair orientations is required, please raise an enchancement request at https://github.com/PapenfussLab/gridss/issues"; 
					    			log.error(msg);
					    			throw new RuntimeException(msg);
					    		}
								input.ensureExtracted();
							} catch (Exception e) {
								log.error(e, "Fatal exception thrown by worker thread.");
								if (getContext().getConfig().terminateOnFirstError) {
									System.exit(1);
								}
								throw e;
							}
							return null;
						}
					};
				}
			}))) {
			// throw exception from worker thread here
			future.get();
		}
		log.info("Evidence extraction complete.");
	}
	private void callVariants(ExecutorService threadpool) throws IOException, InterruptedException, ExecutionException {
		File rawCalls = getContext().getFileSystemContext().getBreakpointVcf(OUTPUT);
		if (!OUTPUT.exists()) {
			if (!rawCalls.exists()) {
				IdentifyVariants iv = new IdentifyVariants();
				copyInputs(iv);
				iv.OUTPUT_VCF = rawCalls;
				execute(iv, threadpool);
			}
			AnnotateVariants annVariants = new AnnotateVariants();
			copyInputs(annVariants);
			annVariants.INPUT_VCF = rawCalls;
			annVariants.OUTPUT_VCF = OUTPUT;
			execute(annVariants, threadpool);
		} else {
			String msg = "Error writing variant calls to " + OUTPUT.getAbsolutePath() + ". File already exists. "
					+ "Please delete OUTPUT file."; 
			log.error(msg);
			throw new IOException(msg);
		}
		if (gridss.Defaults.DELETE_TEMPORARY_FILES) {
			FileHelper.delete(rawCalls, true);
		}
	}
	private void execute(MultipleSamFileCommandLineProgram program, ExecutorService threadpool) throws IOException, InterruptedException, ExecutionException {
		int result = program.doWork(threadpool);
		if (result != 0) throw new RuntimeException("Error executing " + program.getClass().getName() + " return status: " + Integer.toString(result));
	}
	public static void main(String[] argv) {
        System.exit(new CallVariants().instanceMain(argv));
    }
	@Override
	public int doWork(ExecutorService threadpool) throws IOException, InterruptedException, ExecutionException {
		IOUtil.assertFileIsWritable(OUTPUT);
    	extractEvidence(threadpool, getSamEvidenceSources());
    	AssemblyEvidenceSource assemblyEvidence = new AssemblyEvidenceSource(getContext(), getSamEvidenceSources(), ASSEMBLY);
    	if (!ASSEMBLY.exists()) {
    		assemblyEvidence.assembleBreakends(threadpool);
    	}
    	// convert breakend assemblies into breakpoint via split read identification
    	assemblyEvidence.ensureExtracted();
    	// call and annotate variants
    	callVariants(threadpool);
		return 0;
	}
}
