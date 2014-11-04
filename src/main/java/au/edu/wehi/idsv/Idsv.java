package au.edu.wehi.idsv;

import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.EnumSet;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import picard.cmdline.Option;
import picard.cmdline.Usage;
import au.edu.wehi.idsv.pipeline.SortRealignedAssemblies;
import au.edu.wehi.idsv.pipeline.SortRealignedSoftClips;
import au.edu.wehi.idsv.validation.TruthAnnotator;
import au.edu.wehi.idsv.vcf.VcfSvConstants;

import com.google.common.base.Function;
import com.google.common.collect.Lists;
import com.google.common.util.concurrent.MoreExecutors;

/**
 * Extracts structural variation evidence and assembles breakends
 * @author Daniel Cameron
 *
 */
public class Idsv extends CommandLineProgram {
	private Log log = Log.getInstance(Idsv.class);
	private static final String PROGRAM_VERSION = "0.1";
	@Option(doc="Number of worker threads to spawn."
			+ "CPU usage can be higher than the number of worker thread due to additional I/O threads",
    		shortName="THREADS")
    public int WORKER_THREADS = Runtime.getRuntime().availableProcessors();
    // The following attributes define the command-line arguments
	@Option(doc="VCF containing true calls. Called variants will be annotated as true/false positive calls based on this truth file.", optional=true)
    public File TRUTH_VCF = null;
    @Usage
    public String USAGE = getStandardUsagePreamble() + "Calls structural variations from NGS sequencing data. " + PROGRAM_VERSION;
    @Option(doc = "Processing steps to execute",
            optional = true)
	public EnumSet<ProcessStep> STEPS = ProcessStep.ALL_STEPS;
    private SAMEvidenceSource constructSamEvidenceSource(File file, int index, boolean isTumour) throws IOException { 
    	List<Integer> maxFragSizeList = isTumour ? INPUT_TUMOUR_READ_PAIR_MAX_CONCORDANT_FRAGMENT_SIZE : INPUT_READ_PAIR_MAX_CONCORDANT_FRAGMENT_SIZE;
    	if (maxFragSizeList != null && maxFragSizeList.size() > index && maxFragSizeList.get(index) != null) {
    		return new SAMEvidenceSource(getContext(), file, isTumour, 0, maxFragSizeList.get(index));
    	} else if (READ_PAIR_CONCORDANT_PERCENT != null) {
    		return new SAMEvidenceSource(getContext(), file, isTumour, READ_PAIR_CONCORDANT_PERCENT);
    	} else {
    		return new SAMEvidenceSource(getContext(), file, isTumour);
    	}
    }
    @Override
	protected int doWork() {
    	try {
    		ExecutorService threadpool = Executors.newFixedThreadPool(WORKER_THREADS);
    		log.info(String.format("Using %d worker threads", WORKER_THREADS));
	    	ensureDictionariesMatch();

	    	List<SAMEvidenceSource> samEvidence = Lists.newArrayList();
	    	int i = 0;
	    	for (File f : INPUT) {
	    		log.debug("Processing normal input ", f);
	    		SAMEvidenceSource sref = constructSamEvidenceSource(f, i++, false);
	    		samEvidence.add(sref);
	    	}
	    	for (File f : INPUT_TUMOUR) {
	    		log.debug("Processing tumour input ", f);
	    		SAMEvidenceSource sref = constructSamEvidenceSource(f, i++, true);
	    		samEvidence.add(sref);
	    	}
	    	log.debug("Extracting evidence.");
	    	for (Future<Void> future : threadpool.invokeAll(Lists.transform(samEvidence, new Function<SAMEvidenceSource, Callable<Void>>() {
					@Override
					public Callable<Void> apply(final SAMEvidenceSource input) {
						return new Callable<Void>() {
							@Override
							public Void call() throws Exception {
								try {
									input.completeSteps(EnumSet.of(ProcessStep.CALCULATE_METRICS, ProcessStep.EXTRACT_SOFT_CLIPS, ProcessStep.EXTRACT_READ_PAIRS));
								} catch (Exception e) {
									log.error(e, "Exception thrown by worker thread");
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
	    	log.debug("Evidence extraction complete.");
	    	AssemblyEvidenceSource aes = new AssemblyEvidenceSource(getContext(), samEvidence, OUTPUT);
	    	aes.ensureAssembled(threadpool);
	    	
	    	List<EvidenceSource> allEvidence = Lists.newArrayList();
	    	allEvidence.add(aes);
	    	allEvidence.addAll(samEvidence);
	    	
	    	String instructions = getRealignmentScript(allEvidence);
	    	if (instructions != null && instructions.length() > 0) {
	    		log.error("Please realign intermediate fastq files. Suggested command-line for alignment is:\n" +
	    				"##################################\n"+
	    				instructions +
	    				"##################################");
		    	log.error("Please rerun after alignments have been performed.");
		    	if (SCRIPT != null) {
		    		FileWriter writer = null;
		    		try {
		    			writer = new FileWriter(SCRIPT);
		    			writer.write(instructions);
		    		} finally {
		    			if (writer != null) writer.close(); 
		    		}
		    		log.error("Realignment script has been written to ", SCRIPT);
		    	}
		    	return -1;
	    	}
	    	for (Future<Void> future : threadpool.invokeAll(Lists.transform(samEvidence, new Function<SAMEvidenceSource, Callable<Void>>() {
				@Override
				public Callable<Void> apply(final SAMEvidenceSource input) {
					return new Callable<Void>() {
						@Override
						public Void call() throws Exception {
							try {
								SortRealignedSoftClips sortSplitReads = new SortRealignedSoftClips(getContext(), input);
					    		if (!sortSplitReads.isComplete()) {
					    			sortSplitReads.process(STEPS);
					    		}
					    		sortSplitReads.close();
							} catch (Exception e) {
								log.error(e, "Exception thrown by worker thread");
								throw e;
							}
							return null;
						}
					};}}))) {
	    		// throw exception from worker thread here
	    		future.get();
	    	}
	    	SortRealignedAssemblies sra = new SortRealignedAssemblies(getContext(), aes);
    		if (sra.canProcess() && !sra.isComplete()) {
    			sra.process(STEPS, threadpool);
    		}
    		sra.close();
	    	
	    	// check that all steps have been completed
	    	for (SAMEvidenceSource sref : samEvidence) {
	    		if (!sref.isComplete(ProcessStep.CALCULATE_METRICS)
	    			|| !sref.isComplete(ProcessStep.EXTRACT_SOFT_CLIPS)
	    			|| !sref.isComplete(ProcessStep.EXTRACT_READ_PAIRS)
	    			|| !sref.isComplete(ProcessStep.REALIGN_SOFT_CLIPS)
	    			|| !sref.isComplete(ProcessStep.SORT_REALIGNED_SOFT_CLIPS)) {
	    			log.error("Unable to call variants: processing and realignment for ", sref.getSourceFile(), " is not complete.");
	    			return -1;
	    		}
	    	}
	    	if (!aes.isRealignmentComplete()) {
	    		log.error("Unable to call variants: generation and realignment of assemblies not complete.");
    			return -1;
	    	}
	    	
	    	VariantCaller caller = null;
	    	try {
	    		caller = new VariantCaller(getContext(), OUTPUT, allEvidence);
	    		caller.callBreakends(threadpool);
	    		BreakendAnnotator annotator = null;
	    		if (TRUTH_VCF != null) {
	    			annotator = new TruthAnnotator(getContext(), TRUTH_VCF);	
	    		}
	    		caller.annotateBreakpoints(annotator);
	    	} finally {
	    		if (caller != null) caller.close();
	    	}
	    	hackOutputIndels();
    	} catch (IOException e) {
    		log.error(e);
    		throw new RuntimeException(e);
    	} catch (InterruptedException e) {
    		log.error("Interrupted waiting for task to complete", e);
    		throw new RuntimeException(e);
		} catch (ExecutionException e) {
			log.error("Exception thrown from background task", e.getCause());
    		throw new RuntimeException("Exception thrown from background task", e.getCause());
		}
		return 0;
    }
	private void hackOutputIndels() throws IOException {
		log.info("DEBUGGING HACK: naive translation to INDEL calls");
		VariantContextDirectedEvidenceIterator it = new VariantContextDirectedEvidenceIterator(getContext(), null, new VCFFileReader(OUTPUT, false).iterator(), null);
		List<IdsvVariantContext> out = Lists.newArrayList();
		while (it.hasNext()) {
			VariantContextDirectedEvidence e = it.next();
			IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(getContext(), e);
			if (e instanceof VariantContextDirectedBreakpoint) {
				VariantContextDirectedBreakpoint bp = (VariantContextDirectedBreakpoint)e;
				BreakpointSummary loc = bp.getBreakendSummary();
				int indelSize = 0;
				if (loc.referenceIndex == loc.referenceIndex2) {
					if (loc.start < loc.start2 && loc.direction == BreakendDirection.Forward && loc.direction2 == BreakendDirection.Backward) {
						// low position of indel
						indelSize = loc.start - loc.start2;
					} else if (loc.start2 < loc.start && loc.direction2 == BreakendDirection.Forward && loc.direction == BreakendDirection.Backward) {
						indelSize = loc.start2 - loc.start;
					}
					indelSize += bp.getBreakpointSequenceString().length() + 1;
				}
				if (indelSize != 0) {
					builder.attribute(VcfSvConstants.SV_LENGTH_KEY, indelSize);
					builder.attribute(VcfSvConstants.SV_TYPE_KEY, indelSize < 0 ? "DEL" : "INS");
					builder.attribute("ALT", bp.getAlternateAllele(0).getDisplayString());
					builder.alleles("N", indelSize < 0 ? "<DEL>" : "<INS>");
					builder.start(Math.min(loc.start, loc.start2));
					builder.stop(Math.max(loc.start, loc.start2));
					builder.attribute(VCFConstants.END_KEY, Math.max(loc.start, loc.start2));
				}
			}
			out.add(builder.make());
		}
		Collections.sort(out, IdsvVariantContext.ByLocationStart);
		VariantContextWriter writer = getContext().getVariantContextWriter(new File(OUTPUT.toString() + ".indel.vcf"));
		for (VariantContext vc : out) {
			writer.add(vc);
		}
		writer.close();
	}
	public static void main(String[] argv) {
        System.exit(new Idsv().instanceMain(argv));
    }
}
