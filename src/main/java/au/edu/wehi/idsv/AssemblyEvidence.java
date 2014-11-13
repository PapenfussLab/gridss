package au.edu.wehi.idsv;

import java.util.List;

import au.edu.wehi.idsv.vcf.VcfFilter;

/**
 * An assembly of underlying reads supporting a breakend
 * @author Daniel Cameron
 *
 */
public interface AssemblyEvidence extends DirectedEvidence {
	/**
	 * Aequence of assembly as if mapped to positive strand of the local breakend location 
	 * @return
	 */
	byte[] getAssemblySequence();
	byte[] getAssemblyAnchorSequence();
	int getAssemblyAnchorLength();
	int getAssemblyBaseCount(EvidenceSubset subset);
	int getAssemblySupportCountReadPair(EvidenceSubset subset);
	int getAssemblyReadPairLengthMax(EvidenceSubset subset);
	int getAssemblySupportCountSoftClip(EvidenceSubset subset);
	int getAssemblySoftClipLengthTotal(EvidenceSubset subset);
	int getAssemblySoftClipLengthMax(EvidenceSubset subset);
	boolean isAssemblyFiltered();
	List<String> getAssemblyFilters();
	void filterAssembly(VcfFilter reason);
}
