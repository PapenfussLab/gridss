package gridss.cmdline;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Locale;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

import org.apache.commons.configuration.ConfigurationException;

import com.google.common.collect.Lists;
import com.google.common.util.concurrent.ThreadFactoryBuilder;

import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.SAMEvidenceSource;
import au.edu.wehi.idsv.alignment.AlignerFactory;
import au.edu.wehi.idsv.configuration.GridssConfiguration;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SequenceUtil;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.sam.CreateSequenceDictionary;

public abstract class MultipleSamFileCommandLineProgram extends ReferenceCommandLineProgram {
	private static final Log log = Log.getInstance(MultipleSamFileCommandLineProgram.class);
	@Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Coordinate-sorted input BAM file.")
    public List<File> INPUT;
	@Option(shortName="IN", doc="Name-sorted input BAM file. This is required for if multiple alignment are reported for each read.", optional=true)
    public List<File> INPUT_NAME_SORTED;
	@Option(doc="Input label. Variant calling evidence breakdowns are reported for each label."
			+ " Default labels correspond to INPUT filenames. "
			+ "When specifying labels, labels must be provided for all input files.", optional=true)
    public List<String> INPUT_LABEL;
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
    		+ "performance. For human WGS, the ENCODE DAC blacklist is recommended.", optional=true)
    public File BLACKLIST = null;
    // --- evidence filtering parameters ---
    @Option(shortName="C", doc = "gridss configuration file containing overrides", optional=true)
    public File CONFIGURATION_FILE = null;
	@Option(doc="Number of worker threads to spawn. Defaults to number of cores available."
			+ " Note that I/O threads are not included in this worker thread count so CPU usage can be higher than the number of worker thread.",
    		shortName="THREADS")
    public int WORKER_THREADS = Runtime.getRuntime().availableProcessors();
	
	private List<SAMEvidenceSource> samEvidence = null;
    private SAMEvidenceSource constructSamEvidenceSource(File file, File nameSortedFile, String label, int minFragSize, int maxFragSize) {
    	int category = getContext().registerCategory(label);
    	if (maxFragSize > 0) {
    		return new SAMEvidenceSource(getContext(), file, nameSortedFile, category, minFragSize, maxFragSize);
    	} else if (READ_PAIR_CONCORDANT_PERCENT != null) {
    		return new SAMEvidenceSource(getContext(), file, nameSortedFile, category, READ_PAIR_CONCORDANT_PERCENT);
    	} else {
    		return new SAMEvidenceSource(getContext(), file, nameSortedFile, category);
    	}
    }
    private <T> T getOffset(List<T> list, int offset, T defaultValue) {
    	if (offset >= list.size()) return defaultValue;
    	T obj = list.get(offset);
    	if (obj == null) return defaultValue;
    	return obj;
    }
    public List<SAMEvidenceSource> getSamEvidenceSources() {
    	if (samEvidence == null) {
	    	if (INPUT_LABEL == null || INPUT_LABEL.size() == 0) {
	    		INPUT_LABEL = INPUT.stream().map(x -> x.getName()).collect(Collectors.toList());
			}
	    	samEvidence = Lists.newArrayList();
	    	for (int i = 0; i < INPUT.size(); i++) {
	    		samEvidence.add(constructSamEvidenceSource(
	    				getOffset(INPUT, i, null),
	    				getOffset(INPUT_NAME_SORTED, i, null),
	    				getOffset(INPUT_LABEL, i, ""),
	    				getOffset(INPUT_MIN_FRAGMENT_SIZE, i, 0),
	    				getOffset(INPUT_MAX_FRAGMENT_SIZE, i, 0)));
	    	}
    	}
    	return samEvidence;
    }
    public void setSamEvidenceSources(List<SAMEvidenceSource> sources) {
    	samEvidence = sources;
    }
    private void ensureArgs() {
		IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
		for (File f : INPUT) {
			IOUtil.assertFileIsReadable(f);
		}
		if (WORKING_DIR == null) {
			WORKING_DIR = new File(".");
		}
		IOUtil.assertDirectoryIsWritable(WORKING_DIR);
	}
    public static void ensureIndexed(File fa) throws IOException {
    	try (ReferenceSequenceFile reference = ReferenceSequenceFileFactory.getReferenceSequenceFile(fa)) {
    		if (!reference.isIndexed()) {
    			String msg = String.format("Unable to find index for %1$s. Please run 'samtools faidx %1$s' or picard tools BuildBamIndex to generate an index file.", fa);
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
			ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(REFERENCE_SEQUENCE);
			SAMSequenceDictionary dictionary = ref.getSequenceDictionary();
			final SamReaderFactory samFactory = SamReaderFactory.makeDefault();
			for (File f : INPUT) {
				SamReader reader = null;
				try {
					reader = samFactory.open(f);
					SequenceUtil.assertSequenceDictionariesEqual(reader.getFileHeader().getSequenceDictionary(), dictionary, f, REFERENCE_SEQUENCE);
				} catch (htsjdk.samtools.util.SequenceUtil.SequenceListsDifferException e) {
					log.error("Reference genome used by ", f, " does not match reference genome ", REFERENCE_SEQUENCE, ". ",
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
    @Override
	protected int doWork() {
    	log.debug("htsjdk buffer size = ", htsjdk.samtools.Defaults.BUFFER_SIZE);
    	log.debug("htsjdk async sam read = ", htsjdk.samtools.Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS);
    	log.debug("htsjdk async sam write = ", htsjdk.samtools.Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS);
    	log.debug("Setting language-neutral locale");
    	java.util.Locale.setDefault(Locale.ROOT);
		// Spam output with gridss parameters used
		getContext();
		// *.sv.bam needs to be indexed if we are to do multi-threaded processing 
		SAMFileWriterFactory.setDefaultCreateIndexWhileWriting(true);
		// Force loading of aligner up-front
		log.info("Loading aligner");
		AlignerFactory.create();
		// TODO: how do we make make sure our log message is printed before we crash the JVM
		// with an IllegalInstruction when using libsswjni on an old computer? 
		//aligner.align_smith_waterman("ACGT".getBytes(), "ACGT".getBytes());
		
    	ensureArgs();
    	ExecutorService threadpool = null;
    	try {
    		ensureIndexed(REFERENCE_SEQUENCE);
    		ensureSequenceDictionary(REFERENCE_SEQUENCE);
    		ensureDictionariesMatch();
    		log.info(String.format("Using %d worker threads", WORKER_THREADS));
    		threadpool = Executors.newFixedThreadPool(getContext().getWorkerThreadCount(), new ThreadFactoryBuilder().setDaemon(false).setNameFormat("Worker-%d").build());
    		return doWork(threadpool);
		} catch (IOException | InterruptedException | ExecutionException e) {
			log.error(e);
    		throw new RuntimeException(e);
		} finally {
			shutdownPool(threadpool);
		}
    }
    public abstract int doWork(ExecutorService threadpool) throws IOException, InterruptedException, ExecutionException;
	private ProcessingContext processContext = null;
	public ProcessingContext getContext() {
		if (processContext == null) {
			GridssConfiguration config;
			try {
				config = new GridssConfiguration(CONFIGURATION_FILE, WORKING_DIR);
			} catch (ConfigurationException e) {
				throw new RuntimeException(e);
			}
			processContext = new ProcessingContext(getFileSystemContext(), REFERENCE_SEQUENCE, null, getDefaultHeaders(), config);
			processContext.setCommandLineProgram(this);
			processContext.setFilterDuplicates(IGNORE_DUPLICATES);
			processContext.setWorkerThreadCount(WORKER_THREADS);
			if (BLACKLIST != null) {
				try {
					processContext.setBlacklist(BLACKLIST);
				} catch (IOException e) {
					log.error(e, "Error loading BED blacklist. ", BLACKLIST);
					throw new RuntimeException(e);
				}
			}
		}
		return processContext;
	}
	public void setContext(ProcessingContext processContext) {
		this.processContext = processContext;
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
	@Override
	protected String[] customCommandLineValidation() {
		String[] val = multiSamFileInputCustomCommandLineValidation();
		if (val != null) return val;
		return super.customCommandLineValidation();
	}
	public String[] multiSamFileInputCustomCommandLineValidation() {
		if (WORKER_THREADS < 1) {
			return new String[] { "WORKER_THREADS must be at least one." };
		}
		if (INPUT == null || INPUT.size() == 0) {
			return new String[] { "No INPUT files specified." };
    	}
		if (INPUT_NAME_SORTED != null && INPUT_NAME_SORTED.size() > 0 && INPUT_NAME_SORTED.size() != INPUT.size()) {
    		return new String[] { "INPUT_NAME_SORTED must omitted or specified for every INPUT." };
    	}
    	if (INPUT_LABEL != null && INPUT_LABEL.size() > 0 && INPUT_LABEL.size() != INPUT.size()) {
    		return new String[] { "INPUT_CATEGORY must omitted or specified for every INPUT." };
    	}
    	if (INPUT_LABEL != null && INPUT_LABEL.stream().anyMatch(x -> x == null || x.equals(""))) {
    		return new String[] { "INPUT_CATEGORY must omitted or specified for every INPUT." };
    	}
    	return null;
	}
	@Override
	public void copyInputs(CommandLineProgram cmd) {
		super.copyInputs(cmd);
		if (cmd instanceof MultipleSamFileCommandLineProgram) {
			MultipleSamFileCommandLineProgram prog = (MultipleSamFileCommandLineProgram) cmd;
			prog.BLACKLIST = BLACKLIST;
			prog.CONFIGURATION_FILE = CONFIGURATION_FILE;
			prog.INPUT = INPUT;
			prog.INPUT_LABEL = INPUT_LABEL;
			prog.INPUT_MAX_FRAGMENT_SIZE = INPUT_MAX_FRAGMENT_SIZE;
			prog.INPUT_MIN_FRAGMENT_SIZE = INPUT_MIN_FRAGMENT_SIZE;
			prog.READ_PAIR_CONCORDANT_PERCENT = READ_PAIR_CONCORDANT_PERCENT;
			prog.WORKER_THREADS = WORKER_THREADS;
			prog.processContext = processContext;
			prog.samEvidence = samEvidence;
		}
	}
}
