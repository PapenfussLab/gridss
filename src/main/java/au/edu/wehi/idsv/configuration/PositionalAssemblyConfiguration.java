package au.edu.wehi.idsv.configuration;

import org.apache.commons.configuration.Configuration;

public class PositionalAssemblyConfiguration {
	public static final String CONFIGURATION_PREFIX = "positional";
	public PositionalAssemblyConfiguration(Configuration config) {
		config = config.subset(CONFIGURATION_PREFIX);
		maxPathLengthMultiple = config.getFloat("maxPathLengthMultiple");
	}
	/**
	 * Maximum length of a single path node. Leaves longer that this length will not be collapsed.
	 * 
	 * This limit is required to ensure that the width of the partial graph loaded into memory
	 * is bounded.    
	 */
	public float maxPathLengthMultiple;
	public int maxPathLengthInBases(int readLength) { return (int)(maxPathLengthMultiple * readLength); }
	public PositionalAssemblyConfiguration(float positionalMaxPathLength) {
		this.maxPathLengthMultiple = positionalMaxPathLength;
	}
}