package au.edu.wehi.idsv;

import java.io.File;

public class AssemblyParameters {
	public AssemblyMethod method = AssemblyMethod.DEBRUIJN_SUBGRAPH;
	/**
	 * De Bruijn graph kmer size
	 */
	public int k = 25;
	/**
	 * Maximum of base mismatches for de bruijn kmer paths to be merged   
	 */
	public int maxBaseMismatchForCollapse = 2;
	/**
	 * Only collapse bubble path with a single entry and exit kmer choice
	 */
	public boolean collapseBubblesOnly = true;
	/**
	 * Allow reuse of reference kmers when assembling subsequent
	 * contigs in an assembly iteration
	 */
	public boolean allReferenceKmerReuse = true;
	/**
	 * Maximum number of contigs per assembly iteration
	 */
	public int maxContigsPerAssembly = 1024;
	/**
	 * Subgraph assembly margin in multiples of max fragment size
	 * 
	 * This determines how long to wait before assembling a subgraph 
	 * Too short and we will assemble a subgraph before all evidence has been added to it
	 * Too long and we will have a greater misassembly rate in repetitive regions
	 */
	public float subgraphAssemblyMargin = 2.5f;
	/**
	 * Maximum number of branches consider at kmer branches.
	 * A value of 1 indicates a greedy traversal
	 */
	public int subgraphAssemblyTraversalMaximumBranchingFactor = Integer.MAX_VALUE;
	/**
	 * Maximum nodes when finding the best path
	 */
	public int maxPathTraversalNodes = Defaults.BEST_PATH_MAX_TRAVERSAL;
	/**
	 * Maximum number of nodes to consider when collapsing paths and leaves
	 */
	//public int maxCollapseTraversalNodes = Defaults.COLLAPSE_PATH_MAX_TRAVERSAL; // TODO: refactor to pass this in
	/**
	 * Directory to save debruijn graph visualisation information to.
	 */
	public File debruijnGraphVisualisationDirectory = null;
	/**
	 * Directory to save debruijn graph visualisation information to.
	 */
	public boolean visualiseAll = Defaults.VISUALISE_ALL;
	/**
	 * Directory to save debruijn graph visualisation information to.
	 */
	public boolean visualiseTimeouts = Defaults.VISUALISE_TIMEOUTS;
	/**
	 * Maximum width (in multiples of maximum fragment size) of subgraph
	 */
	public float maxSubgraphFragmentWidth = Defaults.MAX_SUBGRAPH_WIDTH_IN_FRAGMENT_SIZE_MULTIPLES;
	/**
	 * Minimum number of reads contributing the the assembly
	 */
	public int minReads = 3;
}
