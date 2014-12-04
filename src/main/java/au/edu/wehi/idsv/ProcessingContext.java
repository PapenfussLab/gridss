package au.edu.wehi.idsv;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.fastq.FastqWriterFactory;
import htsjdk.samtools.filter.DuplicateReadFilter;
import htsjdk.samtools.filter.FilteringIterator;
import htsjdk.samtools.metrics.Header;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.MappedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;

import java.io.Closeable;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;

import au.edu.wehi.idsv.util.AutoClosingIterator;
import au.edu.wehi.idsv.vcf.VcfConstants;

import com.google.common.collect.ImmutableList;

/**
 * Processing context for the given record
 * @author Daniel Cameron
 *
 */
public class ProcessingContext implements Closeable {
	private final ReferenceSequenceFile reference;
	private final File referenceFile;
	private final SAMSequenceDictionary dictionary;
	private final LinearGenomicCoordinate linear;
	private final FileSystemContext fsContext;
	private final boolean vcf41mode;
	private final boolean perChr;
	private final AssemblyParameters ap;
	private final SoftClipParameters scp;
	private final ReadPairParameters rrp;
	private final RealignmentParameters rp;
	private final VariantCallingParameters vcp;
	private final List<Header> metricsHeaders;
	private final SAMFileHeader basicHeader;
	private boolean filterDuplicates = true;
	private long calculateMetricsRecordCount = Long.MAX_VALUE; 
	public ProcessingContext(
			FileSystemContext fileSystemContext,
			List<Header> metricsHeaders,
			SoftClipParameters softClipParameters,
			ReadPairParameters readPairParameters,
			AssemblyParameters assemblyParameters,
			RealignmentParameters realignmentParameters,
			VariantCallingParameters variantCallingParameters,
			File ref, boolean perChr, boolean vcf41) {
		this.fsContext = fileSystemContext;
		this.metricsHeaders = metricsHeaders;
		this.scp = softClipParameters;
		this.rrp = readPairParameters;
		this.ap = assemblyParameters;
		this.rp = realignmentParameters;
		this.vcp = variantCallingParameters;
		this.perChr = perChr;
		this.vcf41mode = vcf41;
		this.referenceFile = ref;
		try {
			this.reference = new MappedFastaSequenceFile(this.referenceFile);
		} catch (FileNotFoundException e) {
			throw new RuntimeException("Unabled load fasta " + ref, e);
		}
		if (this.reference.getSequenceDictionary() == null) {
			throw new IllegalArgumentException("Missing sequence dictionary for reference genome. Please create using picard CreateSequenceDictionary.");
		}
		this.dictionary = new DynamicSAMSequenceDictionary(this.reference.getSequenceDictionary());
		this.linear = new LinearGenomicCoordinate(this.dictionary);
		this.basicHeader = new SAMFileHeader();
		this.basicHeader.setSequenceDictionary(this.reference.getSequenceDictionary());
	}
	/**
	 * Creates a new metrics file with appropriate headers for this context 
	 * @return MetricsFile
	 */
	public <A extends MetricBase,B extends Comparable<?>> MetricsFile<A,B> createMetricsFile() {
        final MetricsFile<A,B> file = new MetricsFile<A,B>();
        for (final Header h : metricsHeaders) {
            file.addHeader(h);
        }
        return file;
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
		CloseableIterator<SAMRecord> safeIterator = new AutoClosingIterator<SAMRecord>(rawIterator,  ImmutableList.<Closeable>of(reader));
		//if (isUseAsyncIO()) {
		//	finalIterator = new AsyncBufferedIterator<SAMRecord>(filterIterator, READAHEAD_BUFFERS, READAHEAD_BUFFER_SIZE, file != null ? file.getAbsolutePath() : null);
		//}
		return applyCommonSAMRecordFilters(safeIterator);
	}
	public SAMFileWriterFactory getSamFileWriterFactory(boolean sorted) {
		return new SAMFileWriterFactory()
			.setTempDirectory(fsContext.getTemporaryDirectory())
			.setCreateIndex(sorted);
	}
	/**
	 * Applies filters such as duplicate removal that apply to all SAMRecord parsing
	 * @param iterator raw reads
	 * @return iterator with filtered record excluded
	 */
	public CloseableIterator<SAMRecord> applyCommonSAMRecordFilters(CloseableIterator<SAMRecord> iterator) {
		if (filterDuplicates) {
			iterator = new AutoClosingIterator<SAMRecord>(new FilteringIterator(iterator, new DuplicateReadFilter()), ImmutableList.<Closeable>of(iterator)); 
		}
		return iterator;
	}
	public FastqWriterFactory getFastqWriterFactory(){
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
	public ReferenceSequenceFile getReference() {
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
	public boolean shouldProcessPerChromosome() {
		return perChr;
	}
	/**
	 * Determines whether VCF records should be compatible with VCF v4.1
	 */
	public boolean getVcf41Mode() {
		return vcf41mode;
	}
	@Override
	public void close() throws IOException {
		if (reference != null) reference.close();
	}
	public AssemblyParameters getAssemblyParameters() {
		return ap;
	}
	public SoftClipParameters getSoftClipParameters() {
		return scp;
	}
	public RealignmentParameters getRealignmentParameters() {
		return rp;
	}
	public VariantCallingParameters getVariantCallingParameters() {
		return vcp;
	}
	public ReadPairParameters getReadPairParameters() {
		return rrp;
	}
	public boolean isFilterDuplicates() {
		return filterDuplicates;
	}
	public void setFilterDuplicates(boolean filterDuplicates) {
		this.filterDuplicates = filterDuplicates;
	}
	public long getCalculateMetricsRecordCount() {
		return calculateMetricsRecordCount;
	}
	public void setCalculateMetricsRecordCount(long calculateMetricsRecordCount) {
		this.calculateMetricsRecordCount = calculateMetricsRecordCount;
	}
}
