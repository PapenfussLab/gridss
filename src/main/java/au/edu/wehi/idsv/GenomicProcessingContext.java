package au.edu.wehi.idsv;

import java.io.Closeable;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;

import au.edu.wehi.idsv.bed.IntervalBed;
import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.picard.TwoBitBufferedReferenceSequenceFile;
import au.edu.wehi.idsv.util.AutoClosingIterator;
import au.edu.wehi.idsv.vcf.GridssVcfConstants;
import gridss.cmdline.ReferenceCommandLineProgram;
import htsjdk.samtools.*;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.filter.AggregateFilter;
import htsjdk.samtools.filter.DuplicateReadFilter;
import htsjdk.samtools.filter.FailsVendorReadQualityFilter;
import htsjdk.samtools.filter.FilteringSamIterator;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.filter.SecondaryOrSupplementaryFilter;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import picard.cmdline.CommandLineProgram;

public class GenomicProcessingContext implements Closeable {
	private static final Log log = Log.getInstance(GenomicProcessingContext.class);
	/**
	 * Buffer between chromosomes
	 * Must be greater than VariantCallingParameters.breakendMargin
	 * value this huge helps with debugging as the chromosome index and offset are immediately apparent  
	 */
	public static final long LINEAR_COORDINATE_CHROMOSOME_BUFFER = 10000000000L;
	private ReferenceLookup reference;
	private CommandLineProgram program;
	private final File referenceFile;
	private final SAMSequenceDictionary dictionary;
	private final LinearGenomicCoordinate linear;
	private final FileSystemContext fsContext;
	private final SAMFileHeader basicHeader;
	private File blacklistFile;
	private IntervalBed blacklist;
	private boolean filterDuplicates = true;
	private int workerThreads = 1;
	/**
	 * Create a new genomic processing context
	 * @param fileSystemContext file system context
	 * @param referenceFile reference genome
	 * @param reference Existing loaded reference lookup. If null is passed, the reference will be loaded into memory.
	 */
	public GenomicProcessingContext(FileSystemContext fileSystemContext, File referenceFile, ReferenceLookup reference) {
		this.fsContext = fileSystemContext;
		this.referenceFile = referenceFile;
		if (this.referenceFile != null) {
			ReferenceCommandLineProgram.ensureSequenceDictionary(referenceFile, program);
		}
		if (reference == null) {
			this.reference = LoadSynchronizedReference(referenceFile);
			if (Defaults.ASYNC_CACHE_REFERENCE) {
				BackgroundCacheReference(referenceFile);
			}
		} else {
			this.reference = reference;
		}
		if (this.reference.getSequenceDictionary() == null) {
			throw new IllegalArgumentException("Missing sequence dictionary for reference genome. Please create using picard CreateSequenceDictionary.");
		}
		this.dictionary = this.reference.getSequenceDictionary(); // new DynamicSAMSequenceDictionary(this.reference.getSequenceDictionary());
		this.linear = new PaddedLinearGenomicCoordinate(this.dictionary, LINEAR_COORDINATE_CHROMOSOME_BUFFER, true);
		this.basicHeader = new SAMFileHeader();
		this.basicHeader.setSequenceDictionary(this.reference.getSequenceDictionary());
		addGridssVersionPG(this.basicHeader);

		this.blacklist = new IntervalBed(this.linear);
	}
	private void addGridssVersionPG(SAMFileHeader header) {
		int i = 0;
		while (header.getProgramRecord(String.format("gridss%d", i)) != null) {
			i++;
		}
		SAMProgramRecord pg = new SAMProgramRecord(String.format("gridss%d", i));
		pg.setProgramName("gridss");
		pg.setProgramVersion(this.getClass().getPackage().getImplementationVersion());
		// TODO: walk the stack and extract the command line arguments
		header.addProgramRecord(pg);
	}
	/**
	 * Load a reference genome with synchronized access to prevent threading issues
	 * @param referenceFile reference genome fasta
	 * @return reference genome 
	 */
	@SuppressWarnings("resource")
	protected ReferenceLookup LoadSynchronizedReference(File referenceFile) {
		try {
			ReferenceSequenceFile underlying = new IndexedFastaSequenceFile(referenceFile);
			if (referenceFile.length() > Runtime.getRuntime().maxMemory()) {
				log.error("Caching reference fasta in memory would require more than 50% of the memory allocated to the JVM. Allocate more heap memory to the JVM..");
				throw new RuntimeException("Not enough memory to cache reference fasta.");
			}
			return new TwoBitBufferedReferenceSequenceFile(underlying);
		} catch (FileNotFoundException e) {
			throw new RuntimeException("Unabled load fasta " + referenceFile, e);
		}
	}

	protected void BackgroundCacheReference(final File ref) {
		Thread thread = new Thread(() -> {
			try {
				for (SAMSequenceRecord contig : reference.getSequenceDictionary().getSequences()) {
					// hit each contig once to ensure it's loaded into memory
					reference.getBase(contig.getSequenceIndex(), 1);
				}
			} catch (Exception e) {
				log.error(e, "Background caching of reference genome failed.");
			}
	    });
		thread.setDaemon(true);
		thread.setName("LoadReference");
		thread.start();
	}

	public FileSystemContext getFileSystemContext() {
		return fsContext;
	}

	/**
	 * Gets a reader for the given file
	 * @param file SAM/BAM file
	 * @return  htsjdk reader
	 */
	public SamReader getSamReader(File file) {
		return getSamReaderFactory().open(file);
	}

	public SamReaderFactory getSamReaderFactory() {
		SamReaderFactory factory = SamReaderFactory.makeDefault()
				.referenceSequence(getReferenceFile())
				.validationStringency(ValidationStringency.LENIENT);
				//.enable(Option.INCLUDE_SOURCE_IN_RECORDS); // don't need as we're tracking ourselves using EvidenceSource
		return factory;
	}

	public CloseableIterator<SAMRecord> getSamReaderIterator(SamReader reader) {
		return getSamReaderIterator(reader, null);
	}

	public CloseableIterator<SAMRecord> getSamReaderIterator(SamReader reader, SortOrder expectedOrder) {
		return getSamReaderIterator(reader, expectedOrder, null);
	}

	public CloseableIterator<SAMRecord> getSamReaderIterator(File reader) {
		return getSamReaderIterator(reader, null);
	}

	public CloseableIterator<SAMRecord> getSamReaderIterator(File file, SortOrder expectedOrder) {
		return getSamReaderIterator(getSamReader(file), expectedOrder, file);
	}

	private CloseableIterator<SAMRecord> getSamReaderIterator(SamReader reader, SortOrder expectedOrder, File file) {
		SAMRecordIterator rawIterator = reader.iterator();
		if (expectedOrder != null && expectedOrder != SortOrder.unsorted) {
			rawIterator.assertSorted(expectedOrder);
		}
		// wrap so we're happy to close as many times as we want
		CloseableIterator<SAMRecord> safeIterator = new AutoClosingIterator<SAMRecord>(rawIterator,  reader);
		return applyCommonSAMRecordFilters(safeIterator);
	}

	public SAMFileWriterFactory getSamFileWriterFactory(boolean sorted) {
		return new SAMFileWriterFactory()
			.setTempDirectory(fsContext.getTemporaryDirectory())
			.setCreateIndex(sorted); // also covered by -Dcreate_index=true
	}

	/**
	 * Applies filters such as duplicate removal that apply to all SAMRecord parsing
	 * @param iterator raw reads
	 * @return iterator with filtered record excluded
	 */
	public CloseableIterator<SAMRecord> applyCommonSAMRecordFilters(CloseableIterator<SAMRecord> iterator) {
		return applyCommonSAMRecordFilters(iterator, false);
	}

	/**
	 * Applies filters such as duplicate removal that apply to all SAMRecord parsing
	 * @param iterator raw reads
	 * @param singleAlignmentPerRead should secondary and supplementary alignment be filtered out
	 * @return iterator with filtered record excluded
	 */
	public CloseableIterator<SAMRecord> applyCommonSAMRecordFilters(final CloseableIterator<SAMRecord> iterator, final boolean singleAlignmentPerRead) {
		List<SamRecordFilter> filters = Lists.<SamRecordFilter>newArrayList(new FailsVendorReadQualityFilter());
		if (singleAlignmentPerRead) {
			filters.add(new SecondaryOrSupplementaryFilter());
		}
		if (filterDuplicates) {
			filters.add(new DuplicateReadFilter());
		}
		return new AutoClosingIterator<SAMRecord>(new FilteringSamIterator(iterator, new AggregateFilter(filters)), iterator);
	}

	public VariantContextWriterBuilder getVariantContextWriterBuilder(File output, boolean createIndex) {
		VariantContextWriterBuilder builder = new VariantContextWriterBuilder()
			.setOutputFile(output)
			.setReferenceDictionary(getReference().getSequenceDictionary());
		builder.clearOptions();
		if (createIndex) {
			builder.setOption(Options.INDEX_ON_THE_FLY);
		} else {
			builder.clearIndexCreator();
		}
		return builder;
	}
	/**
	 * Gets a VCF file ready to write variants to
	 * A header based on this processing context will have already been written to the returned writer
	 * It is the responsibility of the caller to close the returned @link {@link VariantContextWriter}
	 * @param file file
	 * @return opened output VCF stream
	 */
	public VariantContextWriter getVariantContextWriter(File file, boolean createIndex) {
		return getVariantContextWriter(file, new VCFHeader(), createIndex);
	}
	/**
	 * Gets a VCF file ready to write variants to
	 * A header based on this processing context will have already been written to the returned writer
	 * It is the responsibility of the caller to close the returned @link {@link VariantContextWriter}
	 * @param file file
	 * @return opened output VCF stream
	 */
	public VariantContextWriter getVariantContextWriter(File file, VCFHeader vcfHeader, boolean createIndex) {
		VariantContextWriterBuilder builder = getVariantContextWriterBuilder(file, createIndex);
		VariantContextWriter vcfWriter = builder.build();
		GridssVcfConstants.addHeaders(vcfHeader);
		vcfHeader.setSequenceDictionary(getReference().getSequenceDictionary());
		vcfWriter.writeHeader(vcfHeader);
		return vcfWriter;
	}
	/**
	 * Gets a basic minimal SAM file header matching the reference sequence
	 * @return
	 */
	public SAMFileHeader getBasicSamHeader() {
		return basicHeader.clone();
	}

	public ReferenceLookup getReference() {
		return reference;
	}

	public File getReferenceFile() {
		return referenceFile;
	}

	public SAMSequenceDictionary getDictionary() {
		return dictionary;
	}

	public LinearGenomicCoordinate getLinear() {
		return linear;
	}

	@Override
	public void close() throws IOException {
		log.debug("close() called");
		if (reference != null) reference.close();
	}

	public boolean isFilterDuplicates() {
		return filterDuplicates;
	}

	public void setFilterDuplicates(boolean filterDuplicates) {
		this.filterDuplicates = filterDuplicates;
	}

	public int getWorkerThreadCount() {
		return workerThreads;
	}

	public void setWorkerThreadCount(int workerThreads) {
		this.workerThreads = workerThreads;
	}

	public IntervalBed getBlacklistedRegions() {
		return blacklist;
	}

	public File getBlacklist() {
		return blacklistFile;
	}
	
	public void setBlacklist(File blacklistFile) throws IOException {
		if (!blacklistFile.exists()) {
			throw new IllegalArgumentException(String.format("Missing file %s", blacklistFile));
		}
		this.blacklistFile = blacklistFile;
		this.blacklist = new IntervalBed(getLinear(), blacklistFile);
	}
	
	public CommandLineProgram getCommandLineProgram() {
		return program;
	}
	
	public void setCommandLineProgram(CommandLineProgram program) {
		this.program = program;
	}
}