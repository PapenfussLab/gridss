package au.edu.wehi.idsv.configuration;

import org.apache.commons.configuration.Configuration;

import java.io.File;

public class VisualisationConfiguration {
	public static final String CONFIGURATION_PREFIX = "visualisation";
	public VisualisationConfiguration(Configuration config, File workingDirectory) {
		config = config.subset(CONFIGURATION_PREFIX);
		directory = new File(workingDirectory, config.getString("directory"));
		timeouts = config.getBoolean("timeouts");
		assemblyGraph = config.getBoolean("assemblyGraph");
		assemblyGraphFullSize = config.getBoolean("assemblyGraphFullSize");
		assemblyProgress = config.getBoolean("assemblyProgress");
		assemblyContigMemoization = config.getBoolean("assemblyContigMemoization");
		assemblyTelemetry = config.getBoolean("assemblyTelemetry");
		maxCliqueTelemetry = config.getBoolean("maxCliqueTelemetry");
		evidenceAllocation = config.getBoolean("evidenceAllocation");
		buffers = config.getBoolean("buffers");
		bufferTrackingItervalInSeconds = config.getFloat("bufferTrackingItervalInSeconds");
		
		if (!directory.exists() && (timeouts || assemblyGraph || assemblyGraphFullSize || assemblyProgress || evidenceAllocation || buffers || maxCliqueTelemetry)) {
			directory.mkdir();
		}
	}
	public File directory;
	public boolean timeouts;
	public boolean assemblyGraph;
	/**
	 * Full-size positional de Bruijn assembly graph. This is a very large graph 
	 */
	public boolean assemblyGraphFullSize;
	/**
	 * Output information on assembly progress and buffer sizes 
	 */
	public boolean assemblyProgress;
	/**
	 * Output contig memoization lookup after each generated contig
	 */
	public boolean assemblyContigMemoization;
	public boolean evidenceAllocation;
	public boolean buffers;
	public float bufferTrackingItervalInSeconds;
	public boolean assemblyTelemetry;
	public boolean maxCliqueTelemetry;
}
