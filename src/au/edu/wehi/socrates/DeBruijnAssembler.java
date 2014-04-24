package au.edu.wehi.socrates;

public class DeBruijnAssembler implements ReadEvidenceAssembler {
	private final int k;
	private int currentReferenceIndex = -1;
	private long currentPosition = -1;
	public DeBruijnAssembler(int kmer) {
		this.k = kmer;
	}
	@Override
	public Iterable<AssemblyEvidence> addEvidence(DirectedBreakpoint evidence) {
		// TODO Auto-generated method stub
		return null;
	}
	@Override
	public Iterable<AssemblyEvidence> endOfEvidence() {
		// TODO Auto-generated method stub
		return null;
	}
}
