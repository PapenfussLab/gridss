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
	public int minMapq;
	/**
	 * Minimum entropy of anchored sequence (in bits) (Shannon entropy)
	 */
	public double minAnchorShannonEntropy;
	/**
	 * Maximum per-input coverage after which reads falling within regions exceeding this coverage are excluded from assembly and variant calling
	 */
	public int maxCoverage;
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
	 * Remove the assembly contribution of a multimapping read from all location except the mapping location with the best assembly
	 */
	public boolean multimappingUniqueAssemblyAllocation;
	/**
	 * Remove the variant contribution of a multimapping read/read pair from all location except the mapping location supporting the best variant
	 */
	public boolean multimappingUniqueVariantAllocation;
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
		minMapq = config.getInt("minMapq");
		maxCoverage = config.getInt("maxCoverage");
		minAnchorShannonEntropy = config.getFloat("minAnchorShannonEntropy");
		dovetailMargin = config.getInt("dovetailMargin");
		terminateOnFirstError = config.getBoolean("terminateOnFirstError");
		chunkSize = config.getInt("chunkSize");
		multimappingUniqueAssemblyAllocation = config.getBoolean("multimappingUniqueAssemblyAllocation");
		multimappingUniqueVariantAllocation = config.getBoolean("multimappingUniqueVariantAllocation");
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
