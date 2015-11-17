package au.edu.wehi.idsv.configuration;

import java.io.File;

import org.apache.commons.configuration.Configuration;

public class VisualisationConfiguration {
	public static final String CONFIGURATION_PREFIX = "visualisation";
	public VisualisationConfiguration(Configuration config, File workingDirectory) {
		config = config.subset(CONFIGURATION_PREFIX);
		directory = new File(workingDirectory, config.getString("directory"));
		timeouts = config.getBoolean("timeouts");
		assembly = config.getBoolean("assembly");
		fullSizeAssembly = config.getBoolean("assemblyFullSize");
		assemblyProgress = config.getBoolean("assemblyProgress");
		evidenceAllocation = config.getBoolean("evidenceAllocation");
		buffers = config.getBoolean("buffers");
		bufferTrackingItervalInSeconds = config.getInt("bufferTrackingItervalInSeconds");
		if (!directory.exists() && (timeouts || assembly || fullSizeAssembly || assemblyProgress || evidenceAllocation || buffers)) {
			directory.mkdir();
		}
	}
	public File directory;
	public boolean timeouts;
	public boolean assembly;
	/**
	 * Full-size positional de Bruijn assembly graph. This is a very large graph 
	 */
	public boolean fullSizeAssembly;
	/**
	 * Output information on assembly progress and buffer sizes 
	 */
	public boolean assemblyProgress;
	public boolean evidenceAllocation;
	public boolean buffers;
	public int bufferTrackingItervalInSeconds;
}
