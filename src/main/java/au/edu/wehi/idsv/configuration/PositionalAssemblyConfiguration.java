package au.edu.wehi.idsv.configuration;

import org.apache.commons.configuration.Configuration;

public class PositionalAssemblyConfiguration {
	public static final String CONFIGURATION_PREFIX = "positional";
	public PositionalAssemblyConfiguration(Configuration config) {
		config = config.subset(CONFIGURATION_PREFIX);
		maxPathLengthMultiple = config.getFloat("maxPathLengthMultiple");
		retainWidthMultiple = config.getFloat("retainWidthMultiple");
		flushWidthMultiple = config.getFloat("flushWidthMultiple");
	}
	/**
	 * Maximum length of a single path node. Leaves longer that this length will not be collapsed.
	 * 
	 * This limit is required to ensure that the width of the partial graph loaded into memory
	 * is bounded.    
	 */
	public float maxPathLengthMultiple;
	/**
	 * Maximum width (in multiples of maximum fragment size) of the loaded subgraph before
	 * the frontier start that are not called.
	 * 
	 * This should be greater than @see AssemblyConfiguration.maxExpectedBreakendLengthMultiple
	 */
	public float retainWidthMultiple = 4;
	/**
	 * Maximum width (in multiples of maximum fragment size) over which to call all contigs.
	 * When the distance from the first contig path start exceeds this flushWidthMultiple +
	 * @see retainWidthMultiple, contigs starting earlier than retainWidthMultiple from the
	 * frontier will be called
	 */
	public float flushWidthMultiple = 10;
	public int maxPathLengthInBases(int readLength) { return (int)(maxPathLengthMultiple * readLength); }
	public PositionalAssemblyConfiguration(float positionalMaxPathLength) {
		this.maxPathLengthMultiple = positionalMaxPathLength;
	}
}