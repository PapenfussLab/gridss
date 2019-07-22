package au.edu.wehi.idsv.configuration;

import org.apache.commons.configuration.Configuration;

public class PositionalAssemblyConfiguration {
	public static final String CONFIGURATION_PREFIX = "positional";
	public PositionalAssemblyConfiguration(Configuration config) {
		config = config.subset(CONFIGURATION_PREFIX);
		maxPathLengthMultiple = config.getFloat("maxPathLengthMultiple");
		retainWidthMultiple = config.getFloat("retainWidthMultiple");
		flushWidthMultiple = config.getFloat("flushWidthMultiple");
		maximumNodeDensity = config.getFloat("maximumNodeDensity");
		trimSelfIntersectingReads = config.getBoolean("trimSelfIntersectingReads");
		forceFullMemoizationRecalculationAt = config.getFloat("forceFullMemoizationRecalculationAt");
		safetyModePathCountThreshold = config.getInt("safetyModePathCountThreshold");
		safetyModeContigsToCall = config.getInt("safetyModeContigsToCall");
		if (retainWidthMultiple < 1) {
			throw new IllegalArgumentException("retainWidthMultiple must be at least 1");
		}
		if (flushWidthMultiple < 1) {
			throw new IllegalArgumentException("flushWidthMultiple must be at least 1");
		}
		if (maximumNodeDensity <= 0) {
			throw new IllegalArgumentException("maximumNodeDensity must be positive");
		}
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
	 * This must be greater than @see AssemblyConfiguration.maxExpectedBreakendLengthMultiple
	 */
	public float retainWidthMultiple;
	/**
	 * Maximum width (in multiples of maximum fragment size) over which to call all contigs.
	 * When the distance from the first contig path start exceeds this flushWidthMultiple +
	 * @see retainWidthMultiple, contigs starting earlier than retainWidthMultiple from the
	 * frontier will be called
	 */
	public float flushWidthMultiple;
	/**
	 * Maximum post-compression node density. Assembly will not be performed on regions of the
	 * genome with a density higher than maximumNodeDensity per base pair. 
	 */
	public float maximumNodeDensity;
	public int maxPathLengthInBases(int readLength) { return (int)(maxPathLengthMultiple * readLength); }
	/**
	 * Removes self-intersecting kmers from reads prior to inclusion in the positional de Bruijn graph.
	 * Removing such nodes reduces the misassembly rate and improves runtime performance.
	 */
	public boolean trimSelfIntersectingReads;
	/**
	 * Force full rememoization recalculation whenever more that this portion of the graph is to
	 * be removed.
	 */
	public float forceFullMemoizationRecalculationAt;
	/**
	 * Number of contigs called in safety mode
	 */
	public final int safetyModeContigsToCall;
	/**
	 * Number of memoized paths to enter safety mode
	 */
	public final int safetyModePathCountThreshold;
}