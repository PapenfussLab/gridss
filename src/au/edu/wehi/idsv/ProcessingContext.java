package au.edu.wehi.idsv;

import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.fastq.FastqWriterFactory;
import htsjdk.samtools.filter.DuplicateReadFilter;
import htsjdk.samtools.filter.FilteringIterator;
import htsjdk.samtools.metrics.Header;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.List;

import au.edu.wehi.idsv.vcf.VcfConstants;

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
	private final RealignmentParameters rp;
	private final VariantCallingParameters vcp;
	private final List<Header> metricsHeaders;
	private boolean filterDuplicates = true;
	private boolean useAsyncIO = true;
	public ProcessingContext(
			FileSystemContext fileSystemContext,
			List<Header> metricsHeaders,
			SoftClipParameters softClipParameters,
			AssemblyParameters assemblyParameters,
			RealignmentParameters realignmentParameters,
			VariantCallingParameters variantCallingParameters,
			File ref, boolean perChr, boolean vcf41) {
		this.fsContext = fileSystemContext;
		this.metricsHeaders = metricsHeaders;
		this.scp = softClipParameters;
		this.ap = assemblyParameters;
		this.rp = realignmentParameters;
		this.vcp = variantCallingParameters;
		this.perChr = perChr;
		this.vcf41mode = vcf41;
		this.referenceFile = ref;
		this.reference = ReferenceSequenceFileFactory.getReferenceSequenceFile(this.referenceFile);
		if (this.reference.getSequenceDictionary() == null) {
			throw new IllegalArgumentException("Missing sequence dictionary for reference genome. Please create using picard CreateSequenceDictionary.");
		}
		this.dictionary = new DynamicSAMSequenceDictionary(this.reference.getSequenceDictionary());
		this.linear = new LinearGenomicCoordinate(this.dictionary);
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
	public SamReaderFactory getSamReaderFactory() {
		SamReaderFactory factory = SamReaderFactory.makeDefault()
				.validationStringency(ValidationStringency.LENIENT);
				//.enable(Option.INCLUDE_SOURCE_IN_RECORDS); // don't need as we're tracking ourselves using EvidenceSource
		return factory;
	}
	public SAMFileWriterFactory getSamReaderWriterFactory() {
		return new SAMFileWriterFactory()
			.setTempDirectory(fsContext.getTemporaryDirectory())
			.setUseAsyncIo(isUseAsyncIO())
			.setCreateIndex(true);
	}
	/**
	 * Applies filters such as duplicate removal that apply to all SAMRecord parsing
	 * @param iterator raw reads
	 * @return iterator with filtered record excluded
	 */
	public CloseableIterator<SAMRecord> applyCommonSAMRecordFilters(CloseableIterator<SAMRecord> iterator) {
		if (filterDuplicates) {
			iterator = new FilteringIterator(iterator, new DuplicateReadFilter()); 
		}
		return iterator;
	}
	public FastqWriterFactory getFastqWriterFactory(){
		FastqWriterFactory factory = new FastqWriterFactory();
		factory.setUseAsyncIo(isUseAsyncIO());
		return factory;
	}
	public VariantContextWriterBuilder getVariantContextWriterBuilder(File output) {
		VariantContextWriterBuilder builder = new VariantContextWriterBuilder()
			.setOutputFile(output)
			.setReferenceDictionary(getReference().getSequenceDictionary());
		if (isUseAsyncIO()) {
			builder.setOption(Options.USE_ASYNC_IO);
		}
		return builder;
	}
	/**
	 * Gets a VCF file ready to write variants to
	 * A header based on this processing context will have already been written to the returned writer
	 * It is the responsibility of the caller to close the returned @link {@link VariantContextWriter}
	 * @param output file
	 * @return
	 */
	public VariantContextWriter getVariantContextWriter(File file) {
		VariantContextWriterBuilder builder = getVariantContextWriterBuilder(file);
		VariantContextWriter vcfWriter = builder.build();
		final VCFHeader vcfHeader = new VCFHeader();
		VcfConstants.addHeaders(vcfHeader);
		vcfHeader.setSequenceDictionary(getReference().getSequenceDictionary());
		vcfWriter.writeHeader(vcfHeader);
		return vcfWriter;
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
	public boolean isUseAsyncIO() {
		return useAsyncIO;
	}
	public void setUseAsyncIO(boolean useAsyncIO) {
		this.useAsyncIO = useAsyncIO;
	}
	public boolean isFilterDuplicates() {
		return filterDuplicates;
	}
	public void setFilterDuplicates(boolean filterDuplicates) {
		this.filterDuplicates = filterDuplicates;
	}
}
