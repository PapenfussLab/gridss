package au.edu.wehi.idsv;

import java.io.File;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;

import au.edu.wehi.idsv.configuration.AssemblyConfiguration;
import au.edu.wehi.idsv.configuration.GridssConfiguration;
import au.edu.wehi.idsv.configuration.SoftClipConfiguration;
import au.edu.wehi.idsv.configuration.VariantCallingConfiguration;
import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.vcf.GridssVcfConstants;
import au.edu.wehi.idsv.visualisation.BufferTracker;
import au.edu.wehi.idsv.visualisation.TrackedBuffer;
import htsjdk.samtools.metrics.Header;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;

/**
 * Processing context for the given record
 * @author Daniel Cameron
 *
 */
public class ProcessingContext extends GenomicProcessingContext {
	//private static final Log log = Log.getInstance(ProcessingContext.class);
	private final GridssConfiguration config;
	private final List<Header> metricsHeaders;
	private long calculateMetricsRecordCount = Long.MAX_VALUE; 
	private final List<String> categories = Lists.newArrayList();
	private EvidenceIdentifierGenerator eidgen;
	private BufferTracker bufferTracker = null;
	
	public ProcessingContext(
			FileSystemContext fileSystemContext,  File ref, ReferenceLookup reference, List<Header> metricsHeaders,
			GridssConfiguration config) {
		super(fileSystemContext, ref, reference);
		this.metricsHeaders = metricsHeaders;
		this.config = config;
		if (config.getVisualisation().buffers) {
			bufferTracker = new BufferTracker(new File(config.getVisualisation().directory, "gridss.buffers.csv"), config.getVisualisation().bufferTrackingItervalInSeconds);
			bufferTracker.start();
		}
		this.eidgen = config.hashEvidenceID ? new HashedEvidenceIdentifierGenerator() : new StringEvidenceIdentifierGenerator();
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
	public GridssConfiguration getConfig() {
		return config;
	}
	public AssemblyConfiguration getAssemblyParameters() {
		return getConfig().getAssembly();
	}
	public SoftClipConfiguration getSoftClipParameters() {
		return getConfig().getSoftClip();
	}
	public VariantCallingConfiguration getVariantCallingParameters() {
		return getConfig().getVariantCalling();
	}
	public long getCalculateMetricsRecordCount() {
		return calculateMetricsRecordCount;
	}
	public void setCalculateMetricsRecordCount(long calculateMetricsRecordCount) {
		this.calculateMetricsRecordCount = calculateMetricsRecordCount;
	}
	public void registerBuffer(String context, TrackedBuffer obj) {
		if (bufferTracker != null) {
			bufferTracker.register(context, obj);
		}
	}
	public int registerCategory(String label) {
		int offset = categories.indexOf(label);
		if (offset < 0) {
			categories.add(label);
			offset = categories.size() - 1;
		}
		return offset;
	}
	/**
	 * Number of categories registered  
	 * @return
	 */
	public int getCategoryCount() {
		return categories.size();
	}
	public String getCategoryLabel(int category) {
		return categories.get(category);
	}
	/**
	 * Gets a VCF file ready to write variants to
	 * A header based on this processing context will have already been written to the returned writer
	 * It is the responsibility of the caller to close the returned @link {@link VariantContextWriter}
	 * @param output file
	 * @return opened output VCF stream
	 */
	@Override
	public VariantContextWriter getVariantContextWriter(File file, boolean createIndex) {
		VariantContextWriterBuilder builder = getVariantContextWriterBuilder(file, createIndex);
		VariantContextWriter vcfWriter = builder.build();
		final VCFHeader vcfHeader = new VCFHeader(Collections.emptySet(), categories);
		GridssVcfConstants.addHeaders(vcfHeader);
		vcfHeader.setSequenceDictionary(getReference().getSequenceDictionary());
		vcfWriter.writeHeader(vcfHeader);
		return vcfWriter;
	}
	public EvidenceIdentifierGenerator getEvidenceIDGenerator() {
		return eidgen;
	}
	public void setEvidenceIDGenerator(EvidenceIdentifierGenerator gen) {
		this.eidgen = gen;
	}
}
