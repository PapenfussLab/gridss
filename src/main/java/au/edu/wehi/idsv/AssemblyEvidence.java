package au.edu.wehi.idsv;

import java.util.Collection;
import java.util.List;

import au.edu.wehi.idsv.vcf.VcfFilter;

/**
 * An assembly of underlying reads supporting a breakend
 * @author Daniel Cameron
 *
 */
public interface AssemblyEvidence extends DirectedEvidence {
	/**
	 * Sequence of assembly as if mapped to positive strand of the local breakend location 
	 * @return
	 */
	byte[] getAssemblySequence();
	byte[] getAssemblyAnchorSequence();
	int getAssemblyAnchorLength();
	int getAssemblySupportCountReadPair();
	int getAssemblySupportCountSoftClip();
	int getAssemblySupportCountSoftClipRemote();
	int getAssemblyReadPairLengthMax();
	int getAssemblySoftClipLengthTotal();
	int getAssemblySoftClipLengthMax();
	int getAssemblySupportCountReadPair(int category);
	int getAssemblySupportCountSoftClip(int category);
	int getAssemblySupportCountSoftClipRemote(int category);
	int getAssemblyReadPairLengthMax(int category);
	int getAssemblySoftClipLengthTotal(int category);
	int getAssemblySoftClipLengthMax(int category);
	boolean isAssemblyFiltered();
	List<VcfFilter> getFilters();
	void filterAssembly(VcfFilter reason);
	boolean isPartOfAssembly(DirectedEvidence evidence);
	Collection<DirectedEvidence> getEvidence();
	/**
	 * Gets evidenceIDs of evidence used to construct the assembly
	 * @return
	 */
	Collection<String> getEvidenceIDs();
}
