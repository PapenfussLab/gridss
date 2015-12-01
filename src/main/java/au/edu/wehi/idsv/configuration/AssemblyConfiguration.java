package au.edu.wehi.idsv.configuration;

import org.apache.commons.configuration.Configuration;

import au.edu.wehi.idsv.AssemblyAlgorithm;
import au.edu.wehi.idsv.RealignedRemoteSAMRecordAssemblyEvidence;
import au.edu.wehi.idsv.RemoteEvidence;
import au.edu.wehi.idsv.SAMRecordAssemblyEvidence;
import au.edu.wehi.idsv.vcf.VcfFilter;

public class AssemblyConfiguration {
	public static final String CONFIGURATION_PREFIX = "assembly";
	public AssemblyConfiguration(Configuration config) {
		config = config.subset(CONFIGURATION_PREFIX);
		anchorRealignment = new AnchorRealignmentConfiguration(config);
		subgraph = new SubgraphAssemblyConfiguration(config);
		errorCorrection = new ErrorCorrectionConfiguration(config);
		downsampling = new DownsamplingConfiguration(config);
		positional = new PositionalAssemblyConfiguration(config);
		method = AssemblyAlgorithm.valueOf(config.getString("method"));
		k = config.getInt("k");
		minReads = config.getInt("minReads");
		includeSoftClips = config.getBoolean("includeSoftClips");
		includeAnomalousPairs = config.getBoolean("includeAnomalousPairs");
		includeRemoteSplitReads = config.getBoolean("includeRemoteSplitReads");
		writeFiltered = config.getBoolean("writeFiltered");
		excludeNonSupportingEvidence = config.getBoolean("excludeNonSupportingEvidence");
		anchorLength = config.getInt("anchorLength");
		maxExpectedBreakendLengthMultiple = config.getFloat("maxExpectedBreakendLengthMultiple");
		trackEvidenceID = config.getBoolean("trackEvidenceID");
	}
	public AnchorRealignmentConfiguration anchorRealignment;
	public SubgraphAssemblyConfiguration subgraph;
	public ErrorCorrectionConfiguration errorCorrection;
	public DownsamplingConfiguration downsampling;
	public PositionalAssemblyConfiguration positional;
	/**
	 * Assembly algorithm to use
	 */
	public AssemblyAlgorithm method;
	/**
	 * De Bruijn graph kmer size
	 */
	public int k;
	/**
	 * Minimum number of reads contributing the the assembly
	 */
	public int minReads;
	/**
	 * Include soft clipped reads in assembly
	 */
	public boolean includeSoftClips;
	/**
	 * Include the mate of reads that map to this location and whose mate is not mapped to the expected location
	 */
	public boolean includeAnomalousPairs;
	/**
	 * Include reads with a soft clip that maps to this location
	 */
	public boolean includeRemoteSplitReads;
	/**
	 * Determines whether filtered assemblies are written to intermediate files
	 */
	public boolean writeFiltered;
	public boolean excludeNonSupportingEvidence;
	/**
	 * Default minimum length in bases of reference sequence anchor assembly
	 * 
	 * Note: a breakend assembly longer than this length will cause reference sequence assembly to be at least as long as the breakend 
	 */
	public int anchorLength = 100;
	/**
	 * Maximum expected length of a breakend assembly.
	 * Assemblies larger than this size are extremely likely to be missassemblies
	 * 
	 * Expected max size is 1.0 for single-sided assembly and 2.0 for assembly from both directions 
	 */
	public float maxExpectedBreakendLengthMultiple = 2.5f;
	/**
	 * Retains evidenceID tracking information after evidence rehydration
	 */
	public boolean trackEvidenceID;
	public boolean applyFilters(SAMRecordAssemblyEvidence evidence) {
		SAMRecordAssemblyEvidence localEvidence = evidence;
		if (evidence instanceof RemoteEvidence) {
			// ensures assembly & matching remote always get filtered in/out together
			localEvidence = ((RealignedRemoteSAMRecordAssemblyEvidence)evidence).asLocal();
		}
		boolean filtered = applyBasicFilters(evidence);
		if (localEvidence.getAssemblyAnchorLength() == 0 && localEvidence.getBreakendSequence().length <= localEvidence.getAssemblyReadPairLengthMax()) {
			evidence.filterAssembly(VcfFilter.ASSEMBLY_TOO_SHORT); // assembly length = 1 read
			filtered = true;
		}
		if (localEvidence.getAssemblySupportCountReadPair() == 0 &&
				localEvidence.getAssemblySupportCountSoftClipRemote() == localEvidence.getAssemblySupportCountSoftClip()) { 
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
