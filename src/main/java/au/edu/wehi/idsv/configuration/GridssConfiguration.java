package au.edu.wehi.idsv.configuration;

import java.io.File;
import java.net.URL;

import org.apache.commons.configuration.CompositeConfiguration;
import org.apache.commons.configuration.Configuration;
import org.apache.commons.configuration.ConfigurationException;
import org.apache.commons.configuration.PropertiesConfiguration;

import au.edu.wehi.idsv.AdapterHelper;

/**
 * Configuration settings container for gridss
 * @author Daniel Cameron
 *
 */
public class GridssConfiguration {
	/**
	 * Adapter sequences to ignore
	 */
	public AdapterHelper adapters;
	/**
	 * Minimum MAPQ of read to considered evidence
	 */
	public int minReadMapq;
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
	 * Ignore file timestamps when comparing intermediate and results files.
	 */
	public boolean ignoreFileTimestamps;
	public int async_bufferCount;
	public int async_bufferSize;
	public AssemblyConfiguration getAssembly() {
		return assembly;
	}
	public RealignmentConfiguration getRealignment() {
		return realignment;
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
	private final AssemblyConfiguration assembly;
	private final RealignmentConfiguration realignment;
	private final SoftClipConfiguration softclip;
	private final VisualisationConfiguration visualisation;
	private final VariantCallingConfiguration variantCalling;
	public GridssConfiguration() throws ConfigurationException {
		this((File)null, new File("."));
	}
	public GridssConfiguration(File configuration, File workingDirectory) throws ConfigurationException {
		this(LoadConfig(configuration), workingDirectory);
	}
	public GridssConfiguration(Configuration config, File workingDirectory) {
		assembly = new AssemblyConfiguration(config);
		visualisation = new VisualisationConfiguration(config, workingDirectory);
		realignment = new RealignmentConfiguration(config);
		softclip = new SoftClipConfiguration(config);
		variantCalling = new VariantCallingConfiguration(config);
		adapters = new AdapterHelper(config.getStringArray("adapter"));
		minReadMapq = config.getInt("minReadMapq");
		minAnchorShannonEntropy = config.getFloat("minAnchorShannonEntropy");
		maxCoverage = config.getInt("maxCoverage");
		dovetailMargin = config.getInt("dovetailMargin");
		ignoreFileTimestamps = config.getBoolean("ignoreFileTimestamps");
		async_bufferCount = config.getInt("async.bufferCount");
		async_bufferSize = config.getInt("async.bufferSize");
	}
	private static Configuration LoadConfig(File configuration) throws ConfigurationException {
		CompositeConfiguration config = new CompositeConfiguration();
		if (configuration != null) {
			PropertiesConfiguration configurationOverride = new PropertiesConfiguration(configuration);
			config.addConfiguration(configurationOverride);
		}
		config.addConfiguration(getDefaultPropertiesConfiguration());
		return config;
	}
	private static PropertiesConfiguration defaultConfigFromJar = null;
	private static PropertiesConfiguration getDefaultPropertiesConfiguration() throws ConfigurationException {
		if (defaultConfigFromJar == null) {
			// Load default configuration from jar resources config file
			URL propFileURL = GridssConfiguration.class.getResource("gridss.properties");
			if (propFileURL != null) {
				// maven packaged
				defaultConfigFromJar = new PropertiesConfiguration(propFileURL);
			} else {
				// test cases directly from source
				defaultConfigFromJar = new PropertiesConfiguration(new File("src/main/resources/gridss.properties"));
			}
		}
		return defaultConfigFromJar;
	}
}
