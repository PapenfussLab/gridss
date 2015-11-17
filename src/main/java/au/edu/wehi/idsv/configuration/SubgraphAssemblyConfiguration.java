package au.edu.wehi.idsv.configuration;

import org.apache.commons.configuration.Configuration;

public class SubgraphAssemblyConfiguration {
	public static final String CONFIGURATION_PREFIX = "subgraph";
	public SubgraphAssemblyConfiguration(Configuration contig) {
		contig = contig.subset(CONFIGURATION_PREFIX);
		allReferenceKmerReuse = contig.getBoolean("allReferenceKmerReuse");
		maxSubgraphFragmentWidth = contig.getFloat("maxSubgraphFragmentWidth");
		minSubgraphWidthForTimeout = contig.getInt("minSubgraphWidthForTimeout");
		maxContigsPerAssembly = contig.getInt("maxContigsPerAssembly");
		assemblyMargin = contig.getFloat("assemblyMargin");
		traversalMaximumBranchingFactor = contig.getInt("traversalMaximumBranchingFactor");
		traveralMaximumPathNodes = contig.getInt("traveralMaximumPathNodes");
	}
	/**
	 * Allow reuse of reference kmers when assembling subsequent
	 * contigs in an assembly iteration
	 */
	public boolean allReferenceKmerReuse;
	/**
	 * Maximum width (in multiples of maximum fragment size) of subgraph
	 */
	public float maxSubgraphFragmentWidth;
	/**
	 * Minimum number of bases before immediate assembly is forced
	 */
	public int minSubgraphWidthForTimeout;
	/**
	 * Maximum number of contigs per assembly iteration
	 */
	public int maxContigsPerAssembly;
	/**
	 * Subgraph assembly margin in multiples of max fragment size
	 * 
	 * This determines how long to wait before assembling a subgraph 
	 * Too short and we will assemble a subgraph before all evidence has been added to it
	 * Too long and we will have a greater misassembly rate in repetitive regions
	 * 
	 * Evidence can be up to 1 read length (SC) or 1 fragment size (read pair) from
	 * breakend position. Maximum detectable breakend size is (twice fragment size - kmer size)
	 * = 4 * fragment size for maximal indel breakpoint
	 */
	public float assemblyMargin;
	/**
	 * Maximum number of branches consider at kmer branches.
	 * A value of 1 indicates a greedy traversal
	 */
	public int traversalMaximumBranchingFactor;
	/**
	 * Maximum nodes when finding the best path
	 * (Safety limit to prevent unbounded exponential runtime
	 * when attempt to path collapse highly collected degenerate subgraphs)
	 */
	public int traveralMaximumPathNodes;
}