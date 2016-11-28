package gridss;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import com.google.common.base.Function;
import com.google.common.collect.Lists;
import com.google.common.io.Files;

import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.SAMEvidenceSource;
import au.edu.wehi.idsv.SingleReadEvidence;
import au.edu.wehi.idsv.VariantCaller;
import htsjdk.samtools.SamPairUtil.PairOrientation;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
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
public class CallVariants extends GridssCommandLineProgram {
	private static final Log log = Log.getInstance(CallVariants.class);
	@Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="VCF structural variation calls.")
    public File OUTPUT;
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
	private void callVariants(ExecutorService threadpool, List<SAMEvidenceSource> samEvidence, AssemblyEvidenceSource assemblyEvidence) throws IOException {
		try (VariantCaller caller = new VariantCaller(getContext(), OUTPUT, samEvidence, assemblyEvidence)) {
			//EvidenceToCsv evidenceDump = null;
			//if (getContext().getConfig().getVisualisation().evidenceAllocation) {
			//	evidenceDump = new EvidenceToCsv(new File(getContext().getConfig().getVisualisation().directory, "evidence.csv"));
			//}
			// Run variant caller single-threaded as we can do streaming calls
			caller.callBreakends(threadpool);
			caller.annotateBreakpoints(threadpool);
		}
	}
	private void writeAssemblyBreakends(AssemblyEvidenceSource assemblyEvidence) throws IOException {
		log.info("Writing breakend assembly support.");
		File breakendfa = new File(OUTPUT.getAbsolutePath() + ".breakend.fa.tmp");
		BufferedOutputStream writer = null;
		try {
			writer = new BufferedOutputStream(new FileOutputStream(breakendfa));
			CloseableIterator<DirectedEvidence> it = assemblyEvidence.iterator();
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
			it.close();
			writer.close();
			Files.move(breakendfa, new File(OUTPUT.getAbsolutePath() + ".breakend.fa"));
		} finally {
			CloserUtil.close(writer);
		}
		log.info("Writing breakend assembly support complete.");
	}
	@Override
	protected String[] customCommandLineValidation() {
		if (REFERENCE_SEQUENCE == null) {
            return new String[]{"Must have a non-null REFERENCE_SEQUENCE"};
        }
		if (WORKER_THREADS < 1) {
			return new String[] { "WORKER_THREADS must be at least one." };
		}
		return super.customCommandLineValidation();
	}
	public static void main(String[] argv) {
        System.exit(new CallVariants().instanceMain(argv));
    }
	@Override
	protected int doWork(ExecutorService threadpool) throws IOException, InterruptedException, ExecutionException {
		IOUtil.assertFileIsWritable(OUTPUT);
		List<SAMEvidenceSource> samEvidence = createSamEvidenceSources();
    	extractEvidence(threadpool, samEvidence);

    	File assemblyFile = getContext().getFileSystemContext().getAssembly(OUTPUT);
    	AssemblyEvidenceSource assemblyEvidence = new AssemblyEvidenceSource(getContext(), samEvidence, assemblyFile);
    	if (!assemblyFile.exists()) {
    		assemblyEvidence.assembleBreakends(threadpool);
    	}
    	// convert breakend assemblies into breakpoint via split read identification
    	assemblyEvidence.ensureExtracted();

    	if (!OUTPUT.exists()) {
	    	callVariants(threadpool, samEvidence, assemblyEvidence);
	    	writeAssemblyBreakends(assemblyEvidence);
    	} else {
    		log.info(OUTPUT.toString() + " exists - not recreating.");
    	}
		return 0;
	}
}
