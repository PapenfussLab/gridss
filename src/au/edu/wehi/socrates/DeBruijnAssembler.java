package au.edu.wehi.socrates;

import com.google.common.collect.Lists;

public class DeBruijnAssembler implements ReadEvidenceAssembler {
	private final int k;
	private int currentReferenceIndex = -1;
	private long currentPosition = -1;
	public DeBruijnAssembler(int kmer) {
		this.k = kmer;
	}
	@Override
	public Iterable<DirectedBreakpointAssembly> addEvidence(DirectedEvidence evidence) {
		return Lists.newArrayList();
	}
	@Override
	public Iterable<DirectedBreakpointAssembly> endOfEvidence() {
		return Lists.newArrayList();
	}
}
