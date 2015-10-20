package au.edu.wehi.idsv.configuration;

import java.io.File;

import au.edu.wehi.idsv.AssemblyAlgorithm;
import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.EvidenceSubset;
import au.edu.wehi.idsv.RealignedRemoteSAMRecordAssemblyEvidence;
import au.edu.wehi.idsv.RemoteEvidence;
import au.edu.wehi.idsv.SAMRecordAssemblyEvidence;
import au.edu.wehi.idsv.vcf.VcfFilter;

public class AssemblyConfiguration {
	public DownsamplingConfiguration downsampling = new DownsamplingConfiguration();
	/**
	 * Assembly algorithm to use
	 */
	public AssemblyAlgorithm method = AssemblyAlgorithm.Positional;
	/**
	 * De Bruijn graph kmer size
	 */
	public int k = 25;
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
	public SubgraphAssemblyConfiguration subgraph = new SubgraphAssemblyConfiguration();
	public ErrorCorrectionConfiguration errorCorrection = new ErrorCorrectionConfiguration();
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
	public boolean includeSoftClips = true;
	/**
	 * Include the mate of reads that map to this location and whose mate is not mapped to the expected location
	 */
	public boolean includeAnomalousPairs = true;
	/**
	 * Include reads with a soft clip that maps to this location
	 */
	public boolean includeRemoteSplitReads = true;
	/**
	 * Output internal assembly state information for debugging purposes
	 */
	public boolean trackAlgorithmProgress = Defaults.VISUALISE_ASSEMBLY_PROGRESS || Defaults.VISUALISE_ALL;
	/**
	 * Determines whether filtered assemblies are written to intermediate files
	 */
	public boolean writeFilteredAssemblies = Defaults.WRITE_FILTERED_ASSEMBLIES;
	public boolean excludeNonSupportingEvidence = Defaults.EXCLUDE_ASSEMBLY_NON_SUPPORTING_EVIDENCE;
	/**
	 * Default minimum length in bases of reference sequence anchor assembly
	 * 
	 * Note: a breakend assembly longer than this length will cause reference sequence assembly to be at least as long as the breakend 
	 */
	public int anchorAssemblyLength = 100;
	/**
	 * Maximum expected length of a breakend assembly.
	 * Assemblies larger than this size are extremely likely to be missassemblies
	 * 
	 * Expected max size is 1.0 for single-sided assembly and 2.0 for assembly from both directions 
	 */
	public float maxExpectedBreakendAssemblyLengthInFragmentMultiples = 2.5f;
	public AnchorRealignmentConfiguration anchorRealignment = new AnchorRealignmentConfiguration();
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
