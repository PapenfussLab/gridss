package au.edu.wehi.idsv.configuration;

import java.io.File;
import java.net.URL;

import org.apache.commons.configuration.CompositeConfiguration;
import org.apache.commons.configuration.Configuration;
import org.apache.commons.configuration.ConfigurationException;
import org.apache.commons.configuration.PropertiesConfiguration;

import com.google.common.collect.Lists;

import au.edu.wehi.idsv.AdapterHelper;
import htsjdk.samtools.util.Log;

/**
 * Configuration settings container for gridss
 * @author Daniel Cameron
 *
 */
public class GridssConfiguration {
	private static final Log log = Log.getInstance(GridssConfiguration.class);
	/**
	 * Adapter sequences to ignore
	 */
	public AdapterHelper adapters;
	/**
	 * Minimum MAPQ to considered uniquely mapped
	 */
	public double minMapq;
	/**
	 * Fall-back mapq when no aligner mapping quality is reported  
	 */
	public int fallbackMapq;
	/**
	 * Fall-back base quality score when no Q2 tag is present
	 */
	public int fallbackBaseq;
	/**
	 * Minimum entropy of anchored sequence (in bits) (Shannon entropy)
	 */
	public double minAnchorShannonEntropy;
	/**
	 * Maximum per-input coverage after which reads falling within regions exceeding this coverage are excluded from assembly and variant calling
	 */
	public int maxCoverage;
	/**
	 * Maximum coverage after which reads falling within regions exceeding this coverage are excluded from assembly and variant calling
	 */
	//public int maxCoverage;
	/**
	 * Number of bases reads can overlap in the incorrect orientation and still be considered
	 * concordantly overlapping. Homology between reference and adapter sequences causes alignments
	 * of short fragments to appear to be in the incorrect orientation. 
	 */
	public int dovetailMargin;
	/**
	 * Terminate the JVM whenever the first thread-fatal error is encountered
	 * instead of deferring until all threads have completed.
	 */
	public boolean terminateOnFirstError;
	/**
	 * Number of bases to allocate to each processing task when performing parallel operations
	 * such as assembly.
	 */
	public int chunkSize;
	/**
	 * Penalty (in base pairs) applied to each change of reference sequence for a given chunk.
	 * Zeroo penalty will makes all chunks (except the last) contain the same number of bases,
	 * and a penalty greater than or equal to the chunk size will ensure that each chunk comes
	 * from a single reference sequence.
	 */
	public int chunkSequenceChangePenalty;
	/**
	 * Use the read group sample name as the category label
	 */
	public boolean useReadGroupSampleNameCategoryLabel;
	/**
	 * Use a hashed evidenceID to save space and prevent read names exceeding the 254 character limit imposed by BAM 
	 */
	public boolean hashEvidenceID;
	public AssemblyConfiguration getAssembly() {
		return assembly;
	}
	public SoftClipConfiguration getSoftClip() {
		return softclip;
	}
	public VisualisationConfiguration getVisualisation() {
		return visualisation;
	}
	public VariantCallingConfiguration getVariantCalling() {
		return variantCalling;
	}
	public ScoringConfiguration getScoring() {
		return scoring;
	}
	private final AssemblyConfiguration assembly;
	private final SoftClipConfiguration softclip;
	private final VisualisationConfiguration visualisation;
	private final VariantCallingConfiguration variantCalling;
	private final ScoringConfiguration scoring;
	public GridssConfiguration() throws ConfigurationException {
		this((File)null, new File("."));
	}
	public GridssConfiguration(File configuration, File workingDirectory) throws ConfigurationException {
		this(LoadConfiguration(configuration), workingDirectory);
	}
	public GridssConfiguration(Configuration config, File workingDirectory) {
		assembly = new AssemblyConfiguration(config);
		visualisation = new VisualisationConfiguration(config, workingDirectory);
		softclip = new SoftClipConfiguration(config);
		variantCalling = new VariantCallingConfiguration(config);
		scoring =  new ScoringConfiguration(config);
		adapters = new AdapterHelper(config.getStringArray("adapter"));
		minMapq = config.getDouble("minMapq");
		fallbackMapq = config.getInt("fallbackMapq");
		fallbackBaseq = config.getInt("fallbackBaseq");
		maxCoverage = config.getInt("maxCoverage");
		minAnchorShannonEntropy = config.getFloat("minAnchorShannonEntropy");
		dovetailMargin = config.getInt("dovetailMargin");
		terminateOnFirstError = config.getBoolean("terminateOnFirstError");
		chunkSize = config.getInt("chunkSize");
		chunkSequenceChangePenalty = config.getInt("chunkSequenceChangePenalty");
		useReadGroupSampleNameCategoryLabel = config.getBoolean("useReadGroupSampleNameCategoryLabel");
		hashEvidenceID = config.getBoolean("hashEvidenceID");
	}
	public static Configuration LoadConfiguration(File configuration) throws ConfigurationException {
		CompositeConfiguration config = new CompositeConfiguration();
		if (configuration != null) {
			PropertiesConfiguration configurationOverride = new PropertiesConfiguration(configuration);
			config.addConfiguration(configurationOverride);
		}
		config.addConfiguration(getDefaultPropertiesConfiguration());
		for (String key : Lists.newArrayList(config.getKeys())) {
			for (String value : config.getStringArray(key)) {
				log.info(String.format("%s=%s", key, value));
			}
		}
		return config;
	}
	private static PropertiesConfiguration defaultConfigFromJar = null;
	private static PropertiesConfiguration getDefaultPropertiesConfiguration() throws ConfigurationException {
		if (defaultConfigFromJar == null) {
			// Load default configuration from jar resources config file
			URL propFileURL = GridssConfiguration.class.getResource("/gridss.properties");
			if (propFileURL != null) {
				// maven packaged
				log.debug("Found gridss.properties in jar.");
				defaultConfigFromJar = new PropertiesConfiguration(propFileURL);
			} else {
				File file = new File("src/main/resources/gridss.properties");
				if (file.exists()) {
					log.debug("Using gridss.properties from source.");
					// test cases load directly from source
					defaultConfigFromJar = new PropertiesConfiguration();
				}
			}
			if (defaultConfigFromJar == null) {
				log.error("Unable to load default properties");
			}
		}
		return defaultConfigFromJar;
	}
}
