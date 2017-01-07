package au.edu.wehi.idsv;

import java.io.Closeable;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Path;
import java.util.List;

import com.google.common.collect.Lists;

import au.edu.wehi.idsv.bed.IntervalBed;
import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.picard.TwoBitBufferedReferenceSequenceFile;
import au.edu.wehi.idsv.util.AutoClosingIterator;
import au.edu.wehi.idsv.vcf.VcfConstants;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.fastq.FastqWriterFactory;
import htsjdk.samtools.filter.AggregateFilter;
import htsjdk.samtools.filter.DuplicateReadFilter;
import htsjdk.samtools.filter.FailsVendorReadQualityFilter;
import htsjdk.samtools.filter.FilteringSamIterator;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.filter.SecondaryOrSupplementaryFilter;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;

public class GenomicProcessingContext implements Closeable {
	private static final Log log = Log.getInstance(GenomicProcessingContext.class);
	/**
	 * Buffer between chromosomes
	 * Must be greater than VariantCallingParameters.breakendMargin
	 * value this huge helps with debugging as the chromosome index and offset are immediately apparent  
	 */
	static final long LINEAR_COORDINATE_CHROMOSOME_BUFFER = 10000000000L;
	private ReferenceLookup reference;
	private final File referenceFile;
	private final SAMSequenceDictionary dictionary;
	private final LinearGenomicCoordinate linear;
	private final FileSystemContext fsContext;
	private final SAMFileHeader basicHeader;
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
		this.dictionary = new DynamicSAMSequenceDictionary(this.reference.getSequenceDictionary());
		this.linear = new PaddedLinearGenomicCoordinate(this.dictionary, LINEAR_COORDINATE_CHROMOSOME_BUFFER, true);
		this.basicHeader = new SAMFileHeader();
		this.basicHeader.setSequenceDictionary(this.reference.getSequenceDictionary());
		this.blacklist = new IntervalBed(this.dictionary, this.linear);
	}
	/**
	 * Ensures that a sequence dictionary exists for the given reference
	 * @param referenceFile reference genome fasta
	 */
	protected static void ensureSeqeunceDictionary(File referenceFile) {
		try {
			ReferenceSequenceFile rsf = new FastaSequenceFile(referenceFile, false);
			Path path = referenceFile.toPath().toAbsolutePath();
			if (rsf.getSequenceDictionary() == null) {
				log.info("Attempting to create sequence dictionary for " + referenceFile);
				Path dictPath = path.resolveSibling(path.getFileName().toString() + htsjdk.samtools.util.IOUtil.DICT_FILE_EXTENSION);
				picard.sam.CreateSequenceDictionary csd = new picard.sam.CreateSequenceDictionary();
				csd.OUTPUT = dictPath.toFile();
				csd.REFERENCE = referenceFile;
				csd.instanceMain(new String[0]);
			}
			rsf.close();
		} catch (Exception e) {
			log.error("Sequence dictionary creation failed. Please create using picard CreateSequenceDictionary.", e);
		}
	}
	/**
	 * Load a reference genome with synchronized access to prevent threading issues
	 * @param referenceFile reference genome fasta
	 * @return reference genome 
	 */
	@SuppressWarnings("resource")
	protected static ReferenceLookup LoadSynchronizedReference(File referenceFile) {
		ensureSeqeunceDictionary(referenceFile);
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
			.setTempDirectory(fsContext.getTemporaryDirectory());
			//.setCreateIndex(sorted); // covered by -Dcreate_index=true
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
	 * @param filterSecondaryAlignment should secondary alignment be filtered out
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

	public FastqWriterFactory getFastqWriterFactory() {
		FastqWriterFactory factory = new FastqWriterFactory();
		return factory;
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
	 * @param output file
	 * @return opened output VCF stream
	 */
	public VariantContextWriter getVariantContextWriter(File file, boolean createIndex) {
		VariantContextWriterBuilder builder = getVariantContextWriterBuilder(file, createIndex);
		VariantContextWriter vcfWriter = builder.build();
		final VCFHeader vcfHeader = new VCFHeader();
		VcfConstants.addHeaders(vcfHeader);
		vcfHeader.setSequenceDictionary(getReference().getSequenceDictionary());
		vcfWriter.writeHeader(vcfHeader);
		return vcfWriter;
	}

	/**
	 * Gets a basic minimal SAM file header matching the reference sequence
	 * @return
	 */
	public SAMFileHeader getBasicSamHeader() {
		return basicHeader;
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

	public void setBlacklistedRegions(IntervalBed blacklist) {
		if (blacklist != null) {
			this.blacklist = blacklist;
		}
	}

}