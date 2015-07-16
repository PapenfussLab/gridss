package au.edu.wehi.idsv;

import java.io.File;

import au.edu.wehi.idsv.vcf.VcfFilter;

public class AssemblyParameters {
	/**
	 * Evidence per base to assemble without filtering
	 */
	public double acceptDensityPortion = 0.5;
	/**
	 * Filter evidence to expect this maximum average evidence per base within window 
	 */
	public double targetEvidenceDensity = 10;
	/**
	 * Minimum window size for density calculation 
	 */
	public int minimumDensityWindowSize = 1000;
	/**
	 * Assembly algorithm to use
	 */
	public AssemblyAlgorithm method = AssemblyAlgorithm.Positional;
	/**
	 * Assemble soft clips at the remote realigned position as well as the read mapping location
	 */
	public boolean includeRemoteSoftClips = true;
	/**
	 * De Bruijn graph kmer size
	 */
	public int k = 25;
	/**
	 * Maximum of base mismatches for de bruijn kmer paths to be merged   
	 */
	public int maxBaseMismatchForCollapse = 4;
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
	 * Minimum number of bases before immediate assembly is forced
	 */
	public int minSubgraphWidthForTimeout = 20000; // mammalian mitochondrial size is 15-17kb
	/**
	 * Maximum number of contigs per assembly iteration
	 */
	public int subgraphMaxContigsPerAssembly = 1024;
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
	public float subgraphAssemblyMargin = 4f;
	/**
	 * Maximum number of branches consider at kmer branches.
	 * A value of 1 indicates a greedy traversal
	 */
	public int subgraphAssemblyTraversalMaximumBranchingFactor = Integer.MAX_VALUE;
	/**
	 * Maximum nodes when finding the best path
	 */
	public int subgraphMaxPathTraversalNodes = Defaults.BEST_PATH_MAX_TRAVERSAL;
	/**
	 * Maximum length of path to collapse
	 * Units are multiples of max support width (ie largest max fragment size)
	 */
	public float positionalMaxPathCollapseLength = 2f;
	public int positionalMaxPathCollapseLengthInBases(int readLength) { return (int)(positionalMaxPathCollapseLength * readLength); }
	/**
	 * Maximum length of a single path node. Leaves longer that this length will not be collapsed.
	 * 
	 * This limit is required to ensure that the width of the partial graph loaded into memory
	 * is bounded.    
	 */
	public float positionalMaxPathLength = 1.1f;
	public int positionalMaxPathLengthInBases(int readLength) { return (int)(positionalMaxPathLength * readLength); }
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
	public boolean trackAlgorithmProgress = Defaults.VISUALISE_ASSEMBLY_PROGRESS || Defaults.VISUALISE_ALL;
	/**
	 * Determines whether filtered assemblies are written to intermediate files
	 */
	public boolean writeFilteredAssemblies = Defaults.WRITE_FILTERED_ASSEMBLIES;
	/**
	 * Perform local Smith-Waterman realignment of assemblies
	 */
	public boolean performLocalRealignment = true;
	public int realignmentWindowSize = Defaults.ASSEMBLY_REALIGNMENT_WINDOW_SIZE;
	public boolean excludeNonSupportingEvidence = Defaults.EXCLUDE_ASSEMBLY_NON_SUPPORTING_EVIDENCE;
	/**
	 * Default minimum length in bases of reference sequence anchor assembly
	 * 
	 * Note: a breakend assembly longer than this length will cause reference sequence assembly to be at least as long as the breakend 
	 */
	public int anchorAssemblyLength = 100;
	public boolean applyFilters(SAMRecordAssemblyEvidence evidence) {
		SAMRecordAssemblyEvidence localEvidence = evidence;
		if (evidence instanceof RemoteEvidence) {
			// ensures assembly & matching remote always get filtered in/out together
			localEvidence = ((RealignedRemoteSAMRecordAssemblyEvidence)evidence).asLocal();
		}
		boolean filtered = applyBasicFilters(evidence);
		if (localEvidence.getAssemblyAnchorLength() == 0 && localEvidence.getBreakendSequence().length <= localEvidence.getAssemblyReadPairLengthMax(EvidenceSubset.ALL)) {
			evidence.filterAssembly(VcfFilter.ASSEMBLY_TOO_SHORT); // assembly length = 1 read
			filtered = true;
		}
		if (localEvidence.getAssemblySupportCountRemote(EvidenceSubset.ALL) > 0 &&
				localEvidence.getAssemblySupportCountRemote(EvidenceSubset.ALL) == localEvidence.getAssemblySupportCountSoftClip(EvidenceSubset.ALL) + localEvidence.getAssemblySupportCountReadPair(EvidenceSubset.ALL)) {
			// assembly is entirely made of remote support with no reads mapping to our location at all
			evidence.filterAssembly(VcfFilter.ASSEMBLY_REMOTE);
			filtered = true;
		}
		return filtered;
	}
	/**
	 * Applies filters that do not required the assembly to be annotated or realigned
	 * @param evidence
	 * @return
	 */
	public boolean applyBasicFilters(SAMRecordAssemblyEvidence evidence) {
		SAMRecordAssemblyEvidence localEvidence = evidence;
		if (evidence instanceof RemoteEvidence) {
			// ensures assembly & matching remote always get filtered in/out together
			localEvidence = ((RealignedRemoteSAMRecordAssemblyEvidence)evidence).asLocal();
		}
		boolean filtered = false;
		if (localEvidence.isReferenceAssembly() ||
				localEvidence.getBreakendSummary() == null ||
				localEvidence.getBreakendSequence().length == 0) {
			evidence.filterAssembly(VcfFilter.REFERENCE_ALLELE);
			filtered = true;
		}
		if (localEvidence.getAssemblySupportCount() < minReads) {
			evidence.filterAssembly(VcfFilter.ASSEMBLY_TOO_FEW_READ);
			filtered = true;
		}
		return filtered;
	}
}
