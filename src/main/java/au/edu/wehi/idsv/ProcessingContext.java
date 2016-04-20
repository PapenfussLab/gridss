package au.edu.wehi.idsv;

import htsjdk.samtools.metrics.Header;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;

import java.io.File;
import java.util.List;

import au.edu.wehi.idsv.configuration.AssemblyConfiguration;
import au.edu.wehi.idsv.configuration.GridssConfiguration;
import au.edu.wehi.idsv.configuration.RealignmentConfiguration;
import au.edu.wehi.idsv.configuration.SoftClipConfiguration;
import au.edu.wehi.idsv.configuration.VariantCallingConfiguration;
import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.visualisation.BufferTracker;
import au.edu.wehi.idsv.visualisation.TrackedBuffer;

import com.google.common.collect.Lists;

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
	private AssemblyIdGenerator assemblyIdGenerator = new SequentialIdGenerator("asm");
	private final List<String> categories = Lists.newArrayList((String)null);
	private BufferTracker bufferTracker = null;
	
	public ProcessingContext(
			FileSystemContext fileSystemContext,  File ref, boolean perChr, ReferenceLookup reference,
			List<Header> metricsHeaders, GridssConfiguration config) {
		super(fileSystemContext, ref, perChr, reference);
		this.metricsHeaders = metricsHeaders;
		this.config = config;
		if (config.getVisualisation().buffers) {
			bufferTracker = new BufferTracker(new File(config.getVisualisation().directory, "gridss.buffers.csv"), config.getVisualisation().bufferTrackingItervalInSeconds);
			bufferTracker.start();
		}
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
	public RealignmentConfiguration getRealignmentParameters() {
		return getConfig().getRealignment();
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
	public AssemblyIdGenerator getAssemblyIdGenerator() {
		return assemblyIdGenerator;
	}
	public void setAssemblyIdGenerator(AssemblyIdGenerator generator) {
		this.assemblyIdGenerator = generator;
	}
	public void registerBuffer(String context, TrackedBuffer obj) {
		if (bufferTracker != null) {
			bufferTracker.register(context, obj);
		}
	}
	public void registerCategories(int categoryCount) {
		for (int i = 0; i < categoryCount; i++) {
			registerCategory(categoryCount, "");
		}
	}
	public void registerCategory(int category, String description) {
		if (category < 0) throw new IllegalArgumentException("Category cannot be negative");
		while (categories.size() <= category) categories.add(null);
		categories.set(category, description);
	}
	/**
	 * Number of categories registered  
	 * @return
	 */
	public int getCategoryCount() {
		return categories.size();
	}
}
