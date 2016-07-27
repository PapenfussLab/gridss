package au.edu.wehi.idsv;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.EnumSet;
import java.util.List;
import java.util.Locale;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

import org.apache.commons.configuration.ConfigurationException;

import com.google.common.base.Function;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.io.Files;
import com.google.common.util.concurrent.ThreadFactoryBuilder;

import au.edu.wehi.idsv.alignment.AlignerFactory;
import au.edu.wehi.idsv.bed.IntervalBed;
import au.edu.wehi.idsv.configuration.GridssConfiguration;
import au.edu.wehi.idsv.pipeline.SortRealignedSoftClips;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamPairUtil.PairOrientation;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SequenceUtil;
import picard.analysis.InsertSizeMetrics;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.sam.CreateSequenceDictionary;

/**
 * Extracts structural variation evidence and assembles breakends
 * @author Daniel Cameron
 *
 */
@CommandLineProgramProperties(
        usage = "Calls structural variations from one or more SAM/BAM input files.",  
        usageShort = "Calls structural variations from NGS sequencing data"
)
public class Idsv extends CommandLineProgram {
	private static final Log log = Log.getInstance(Idsv.class);
	@Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Coordinate-sorted input BAM file.")
    public List<File> INPUT;
	@Option(shortName="IC", doc="Input category. Variant calling evidence is reported from category 1 (default) to the maximum category specified. "
			+ "Specify categories when you require a breakdown of support (eg tumour/normal or multi-sample variant calling). ", optional=true)
    public List<Integer> INPUT_CATEGORY;	
    @Option(doc = "Per input maximum concordant fragment size.", optional=true)
    public List<Integer> INPUT_MAX_FRAGMENT_SIZE;
    @Option(doc = "Per input minimum concordant fragment size.", optional=true)
    public List<Integer> INPUT_MIN_FRAGMENT_SIZE;
    @Option(doc = "Percent of read pairs considered concorant (0.0-1.0). "
    		+ "If this is unset, the SAM proper pair flag is used to determine whether a read is discordantly aligned. "
    		+ "Explicit fragment size specification overrides this setting.", optional=true)
    public Float READ_PAIR_CONCORDANT_PERCENT = 0.995f;
    @Option(shortName="BL", doc = "BED blacklist of regions to ignore. Assembly of regions such as high-coverage centromeric repeats is slow, "
    		+ "and if such regions are to be filtered in downstream analysis anyway, blacklisting those region will improve runtime "
    		+ "performance. For human WGS, the ENCODE DAC blacklist is recommended. Specify \"none\" to use no blacklist.", optional=true)
    public File BLACKLIST = null;
	@Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="VCF structural variation calls.")
    public File OUTPUT;
    @Option(shortName=StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc="Reference used for alignment")
    public File REFERENCE;
    // --- intermediate file parameters ---
    @Option(doc = "Save intermediate results into separate files for each chromosome."
    		+ " Increases the number of intermediate files but allows a greater level of parallelisation.", optional = true)
    public boolean PER_CHR = true;
    @Option(doc = "Directory to place intermediate results directories. Default location is the same directory"
    		+ " as the associated input or output file.", optional = true)
    public File WORKING_DIR = null;
    // --- evidence filtering parameters ---
    @Option(shortName="C", doc = "gridss configuration file containing overrides", optional=true)
    public File CONFIGURATION_FILE = null;
	@Option(doc="Number of worker threads to spawn. Defaults to number of cores available."
			+ " Note that I/O threads are not included in this worker thread count so CPU usage can be higher than the number of worker thread.",
    		shortName="THREADS")
    public int WORKER_THREADS = Runtime.getRuntime().availableProcessors();
    //@Option(doc = "Processing steps to execute", optional = true)
    public EnumSet<ProcessStep> STEPS = ProcessStep.ALL_STEPS;
    private SAMEvidenceSource constructSamEvidenceSource(File file, int category, int minFragSize, int maxFragSize) {
    	if (maxFragSize > 0) {
    		return new SAMEvidenceSource(getContext(), file, category, minFragSize, maxFragSize);
    	} else if (READ_PAIR_CONCORDANT_PERCENT != null) {
    		return new SAMEvidenceSource(getContext(), file, category, READ_PAIR_CONCORDANT_PERCENT);
    	} else {
    		return new SAMEvidenceSource(getContext(), file, category);
    	}
    }
    private <T> T getOffset(List<T> list, int offset, T defaultValue) {
    	if (offset >= list.size()) return defaultValue;
    	T obj = list.get(offset);
    	if (obj == null) return defaultValue;
    	return obj;
    }
    public List<SAMEvidenceSource> createSamEvidenceSources() {
    	List<SAMEvidenceSource> samEvidence = Lists.newArrayList();
    	for (int i = 0; i < INPUT.size(); i++) {
    		samEvidence.add(constructSamEvidenceSource(
    				getOffset(INPUT, i, null),
    				getOffset(INPUT_CATEGORY, i, 0),
    				getOffset(INPUT_MIN_FRAGMENT_SIZE, i, 0),
    				getOffset(INPUT_MAX_FRAGMENT_SIZE, i, 0)));
    	}
    	return samEvidence;
    }
    private void ensureArgs() {
		IOUtil.assertFileIsReadable(REFERENCE);
		for (File f : INPUT) {
			IOUtil.assertFileIsReadable(f);
		}
		IOUtil.assertFileIsWritable(OUTPUT);
		if (WORKING_DIR == null) {
			WORKING_DIR = new File(".");
		}
		IOUtil.assertDirectoryIsWritable(WORKING_DIR);
	}
    public static void ensureIndexed(File fa) throws IOException {
    	try (ReferenceSequenceFile reference = ReferenceSequenceFileFactory.getReferenceSequenceFile(fa)) {
    		if (!reference.isIndexed()) {
    			String msg = String.format("Unable to find index for %1$s. Please run 'samtools faidx %1$s' to generate an index file.", fa);
    			log.error(msg);
    			throw new IOException(msg);
    		} else {
    			log.debug(fa, " is indexed.");
    		}
    	}
    }
    public static void ensureSequenceDictionary(File fa) throws IOException {
    	File output = new File(fa.getAbsolutePath() + ".dict");
    	try (ReferenceSequenceFile reference = ReferenceSequenceFileFactory.getReferenceSequenceFile(fa)) {
	    	SAMSequenceDictionary dictionary = reference.getSequenceDictionary();
	    	if (dictionary == null) {
	    		log.info("Attempting to generate missing sequence dictionary for ", fa);
	    		try {
		    		final SAMSequenceDictionary sequences = new CreateSequenceDictionary().makeSequenceDictionary(fa);
		            final SAMFileHeader samHeader = new SAMFileHeader();
		            samHeader.setSequenceDictionary(sequences);
		            try (SAMFileWriter samWriter = new SAMFileWriterFactory().makeSAMWriter(samHeader, false, output)) {
		            }
	    		} catch (Exception e) {
	    			log.error("Missing sequence dictionary for ", fa, " and creation of sequencing dictionary failed.",
	    					"Please run the Picard tools CreateSequenceDictionary utility to create", output);
	    			throw e;
	    		}
	    		log.info("Created sequence dictionary ", output);
	    	} else {
	    		log.debug("Found sequence dictionary for ", fa);
	    	}
    	}
    }
	public void ensureDictionariesMatch() throws IOException {
		ReferenceSequenceFile ref = null;
		try {
			ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(REFERENCE);
			SAMSequenceDictionary dictionary = ref.getSequenceDictionary();
			final SamReaderFactory samFactory = SamReaderFactory.makeDefault();
			for (File f : INPUT) {
				SamReader reader = null;
				try {
					reader = samFactory.open(f);
					SequenceUtil.assertSequenceDictionariesEqual(reader.getFileHeader().getSequenceDictionary(), dictionary, f, REFERENCE);
				} catch (htsjdk.samtools.util.SequenceUtil.SequenceListsDifferException e) {
					log.error("Reference genome used by ", f, " does not match reference genome ", REFERENCE, ". ",
							"The reference supplied must match the reference used for every input.");
					throw e;
				} finally {
					if (reader != null) reader.close();
				}
			}
		} finally {
			if (ref != null) ref.close();
		}
	}
	public static String getRealignmentScript(Iterable<? extends EvidenceSource> it) {
    	StringBuilder sb = new StringBuilder();
    	for (EvidenceSource source : it) {
    		if (!source.isRealignmentComplete(true)) {
    			sb.append(source.getRealignmentScript());
    		}
    	}
    	return sb.toString();
    }
    @Override
	protected int doWork() {
    	log.debug("Setting language-neutral locale");
    	java.util.Locale.setDefault(Locale.ROOT);
    	if (INPUT == null || INPUT.size() == 0) {
    		log.error("No INPUT files specified.");
    		return -1;
    	}
    	if (BLACKLIST != null && BLACKLIST.getName().equals("none")) {
    		BLACKLIST = null;
    	}
    	if (INPUT_CATEGORY != null && INPUT_CATEGORY.size() > 0 && INPUT_CATEGORY.size() != INPUT.size()) {
    		log.error("INPUT_CATEGORY must omitted or specified for every INPUT.");
    		return -1;
    	}
    	if (INPUT_CATEGORY != null && INPUT_CATEGORY.stream().anyMatch(x -> x == null)) {
    		log.error("INPUT_CATEGORY must omitted or specified for every INPUT.");
    		return -1;
    	}
    	if (INPUT_CATEGORY != null && INPUT_CATEGORY.stream().anyMatch(x -> x <= 0)) {
    		log.error("INPUT_CATEGORY must be positive integers: negative or zero categories are not valid.");
    		return -1;
    	}
    	// GRIDSS using zero based categories internally - transform arg
    	if (INPUT_CATEGORY == null || INPUT_CATEGORY.size() == 0) {
    		INPUT_CATEGORY = INPUT.stream().map(x -> 0).collect(Collectors.toList());
		} else {
			INPUT_CATEGORY = INPUT_CATEGORY.stream().map(x -> x != null ? x - 1 : 0).collect(Collectors.toList());
		}
    	ensureArgs();
    	ExecutorService threadpool = null;
    	try {
    		ensureIndexed(REFERENCE);
    		ensureSequenceDictionary(REFERENCE);
    		ensureDictionariesMatch();
    		// Spam output with gridss parameters used
    		getContext();
    		if (INPUT_CATEGORY.stream().mapToInt(x -> x).distinct().count() < INPUT_CATEGORY.stream().mapToInt(x -> x).max().orElse(0) - 8) {
        		log.warn("Missing a large number of INPUT_CATEGORY indicies. Performance is likely to be degraded since unused categories are still computed and stored.");
        	}
    		// Force loading of aligner up-front
    		log.info("Loading aligner");
    		AlignerFactory.create();
    		log.info("Testing for presence of external aligner");
    		
    		//hackSimpleCalls(processContext);
    		threadpool = Executors.newFixedThreadPool(getContext().getWorkerThreadCount(), new ThreadFactoryBuilder().setDaemon(false).setNameFormat("Worker-%d").build());
    		log.info(String.format("Using %d worker threads", WORKER_THREADS));
	    	
	    	List<SAMEvidenceSource> samEvidence = createSamEvidenceSources();
	    	
	    	extractEvidence(threadpool, samEvidence);
	    	
	    	for (SAMEvidenceSource e : samEvidence) {
	    		InsertSizeMetrics ism = e.getMetrics().getInsertSizeMetrics();
	    		if (ism != null && ism.PAIR_ORIENTATION != PairOrientation.FR) {
	    			log.error("gridss currently supports only FR read pair orientation. If usage with other read pair orientations is required, please raise an enchancement request at https://github.com/PapenfussLab/gridss/issues");
	    			return -1;
	    		}
	    	}
	    	
	    	if (!checkRealignment(samEvidence)) {
	    		return -1;
    		}
	    	sortRealignedSoftClips(threadpool, samEvidence);
	    	
	    	AssemblyEvidenceSource assemblyEvidence = new AssemblyEvidenceSource(getContext(), samEvidence, OUTPUT);
	    	assemblyEvidence.ensureAssembled(threadpool);
	    	
	    	if (!checkRealignment(ImmutableList.of(assemblyEvidence))) {
	    		return -1;
    		}
	    	// edge case: need to ensure assembled again since realignment could have completed
	    	// without invoking external aligner (0 records to realign). 
	    	assemblyEvidence.ensureAssembled(threadpool);
	    	
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
	    	if (!assemblyEvidence.isRealignmentComplete(true)) {
	    		log.error("Unable to call variants: generation of breakend alignment of assemblies not complete.");
    			return -1;
	    	}
			
	    	if (!OUTPUT.exists()) {
		    	callVariants(threadpool, samEvidence, assemblyEvidence);
		    	writeAssemblyBreakends(assemblyEvidence);
	    	} else {
	    		log.info(OUTPUT.toString() + " exists not recreating.");
	    	}
	    	shutdownPool(threadpool);
			threadpool = null;
	    	//hackSimpleCalls(processContext);
	    	processContext.close();
    	} catch (IOException e) {
    		log.error(e);
    		throw new RuntimeException(e);
    	} catch (InterruptedException e) {
    		log.error("Interrupted waiting for task to complete", e);
    		throw new RuntimeException(e);
		} catch (ExecutionException e) {
			log.error("Exception thrown from background task", e.getCause());
    		throw new RuntimeException("Exception thrown from background task", e.getCause());
		} catch (ConfigurationException e) {
			log.error(e);
    		throw new RuntimeException(e);
		} finally {
			shutdownPool(threadpool);
			threadpool = null;
		}
		return 0;
    }
    private ProcessingContext processContext = null;
	public ProcessingContext getContext() {
		if (processContext == null) {
			FileSystemContext fsc = new FileSystemContext(TMP_DIR.get(0), WORKING_DIR, MAX_RECORDS_IN_RAM);
			GridssConfiguration config;
			try {
				config = new GridssConfiguration(CONFIGURATION_FILE, fsc.getIntermediateDirectory(OUTPUT));
			} catch (ConfigurationException e) {
				throw new RuntimeException(e);
			}
			processContext = new ProcessingContext(fsc, REFERENCE, PER_CHR, null, getDefaultHeaders(), config);
			INPUT_CATEGORY.stream().forEach(x -> processContext.registerCategory(x, ""));
			processContext.setWorkerThreadCount(WORKER_THREADS);
			if (BLACKLIST != null) {
				try {
					processContext.setBlacklistedRegions(new IntervalBed(processContext.getDictionary(), processContext.getLinear(), BLACKLIST));
				} catch (IOException e) {
					log.error(e, "Error loading BED blacklist. ", BLACKLIST);
					throw new RuntimeException(e);
				}
			}
		}
		return processContext;
	}
	private void sortRealignedSoftClips(ExecutorService threadpool, List<SAMEvidenceSource> samEvidence) throws InterruptedException, ExecutionException {
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
		log.info("Evidence extraction complete.");
	}
	private void callVariants(ExecutorService threadpool, List<SAMEvidenceSource> samEvidence, AssemblyEvidenceSource assemblyEvidence) throws IOException, ConfigurationException {
		VariantCaller caller = null;
		try {
			EvidenceToCsv evidenceDump = null;
			if (getContext().getConfig().getVisualisation().evidenceAllocation) {
				evidenceDump = new EvidenceToCsv(new File(getContext().getConfig().getVisualisation().directory, "evidence.csv"));
			}
			// Run variant caller single-threaded as we can do streaming calls
			caller = new VariantCaller(
				getContext(),
				OUTPUT,
				samEvidence,
				assemblyEvidence,
				evidenceDump);
			caller.callBreakends(threadpool);
			caller.annotateBreakpoints(threadpool);
		} finally {
			if (caller != null) caller.close();
		}
	}
	private void writeAssemblyBreakends(AssemblyEvidenceSource assemblyEvidence) throws IOException {
		log.info("Writing breakend assembly support.");
		File output = new File(OUTPUT.getAbsolutePath() + ".breakend.fa.tmp");
		BufferedOutputStream writer = null;
		try {
			writer = new BufferedOutputStream(new FileOutputStream(output));
			CloseableIterator<SAMRecordAssemblyEvidence> it = assemblyEvidence.iterator(false, false);
			while (it.hasNext()) {
				SAMRecordAssemblyEvidence ass = it.next();
				writer.write('>');
				writer.write(ass.getEvidenceID().getBytes(StandardCharsets.US_ASCII));
				writer.write('\n');
				writer.write(ass.getAssemblySequence());
				writer.write('\n');
			}
			it.close();
			writer.close();
			Files.move(output, new File(OUTPUT.getAbsolutePath() + ".breakend.fa"));
		} finally {
			CloserUtil.close(writer);
		}
		log.info("Writing breakend assembly support complete.");
	}
//	private void hackSimpleCalls() throws IOException, ConfigurationException {
//		if (OUTPUT.exists()) {
//			File simpleHack = new File(OUTPUT.getParent(), OUTPUT.getName() + ".simple.vcf");
//			if (!simpleHack.exists()) {
//				new BreakendToSimpleCall(getContext()).convert(OUTPUT, simpleHack);
//			}
//		}
//	}
    private boolean checkRealignment(List<? extends EvidenceSource> evidence) throws IOException {
    	String instructions = getRealignmentScript(evidence);
    	if (instructions != null && instructions.length() > 0) {
    		log.error("Please realign intermediate fastq files. Suggested command-line for alignment is:\n" +
    				"##################################\n"+
    				instructions +
    				"##################################");
	    	log.error("Please rerun after alignments have been performed.");
	    	return false;
    	}
    	return true;
	}
	private void shutdownPool(ExecutorService threadpool) {
    	if (threadpool != null) {
			log.debug("Shutting down thread pool.");
			threadpool.shutdownNow();
			log.debug("Waiting for thread pool tasks to complete");
			try {
				if (!threadpool.awaitTermination(10, TimeUnit.MINUTES)) {
					log.error("Tasks did not respond to termination request in a timely manner - outstanding tasks will be forcibly terminated without cleanup.");
				}
			} catch (InterruptedException e) {
				throw new RuntimeException(e);
			}
		}
    }
	public static void main(String[] argv) {
        System.exit(new Idsv().instanceMain(argv));
    }
}
