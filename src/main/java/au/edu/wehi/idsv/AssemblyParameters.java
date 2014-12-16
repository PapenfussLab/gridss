package au.edu.wehi.idsv;

import htsjdk.samtools.util.Log;

import java.io.File;

import au.edu.wehi.idsv.vcf.VcfFilter;

public class AssemblyParameters {
	private static final Log log = Log.getInstance(AssemblyParameters.class);
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
	/**
	 * Include soft clipped reads in assembly
	 */
	public boolean assemble_soft_clips = true;
	/**
	 * Include the mate of reads that map to this location and whose mate is not mapped to the expected location
	 */
	public boolean assemble_read_pairs = true;
	/**
	 * Include reads with a soft clip that maps to this location
	 */
	public boolean assemble_remote_soft_clips = true;
	/**
	 * Output internal assembly state information for debugging purposes
	 */
	public boolean trackAlgorithmProgress = Defaults.VISUALISE_ASSEMBLY_PROGRESS;
	/**
	 * Determines whether filtered assemblies are written to intermediate files
	 */
	public boolean writeFilteredAssemblies = Defaults.WRITE_FILTERED_ASSEMBLIES;
	public boolean applyFilters(AssemblyEvidence evidence) {
		boolean filtered = false;
		if (evidence.getBreakendSequence() == null) {
			log.error("Breakpoint sequence missing for assembly " + evidence.toString());
		}
		if (evidence.getBreakendSequence().length == 0) {
			evidence.filterAssembly(VcfFilter.ASSEMBLY_REF);
			filtered = true;
		}
		if (evidence.getAssemblySupportCountReadPair(EvidenceSubset.ALL) + evidence.getAssemblySupportCountSoftClip(EvidenceSubset.ALL) < minReads) {
			evidence.filterAssembly(VcfFilter.ASSEMBLY_TOO_FEW_READ);
			filtered = true;
		}
		if (evidence.getAssemblyAnchorLength() == 0 && evidence.getBreakendSequence().length <= evidence.getAssemblyReadPairLengthMax(EvidenceSubset.ALL)) {
			// just assembled a single read - not very exciting
			evidence.filterAssembly(VcfFilter.ASSEMBLY_TOO_SHORT); // assembly length = 1 read
			filtered = true;
		}
		return filtered;
	}
}
