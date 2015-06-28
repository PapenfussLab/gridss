package au.edu.wehi.idsv.debruijn.positional;

/**
 * Evidence contribution to the given kmer position of a single piece of evidence
 * @author cameron.d
 *
 */
public class KmerSupportNode implements KmerNode {
	private final int offset;
	private final KmerEvidence evidence;
	public long lastKmer() { return evidence.kmer(offset); }
	public int lastStart() { return evidence.startPosition() + offset; }
	public int lastEnd() { return evidence.endPosition() + offset; }
	public int weight() { return evidence.weight(offset); }
	public boolean isReference() { return evidence.isAnchored(offset); }
	public KmerEvidence evidence() { return evidence; }
	public KmerSupportNode(KmerEvidence evidence, int offset) {
		this.evidence = evidence;
		this.offset = offset;
	}
	@Override
	public String toString() {
		return String.format("[%d-%d]%s %s", lastStart(), lastEnd(), isReference() ? "R" : " ", weight());
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result
				+ ((evidence == null) ? 0 : evidence.hashCode());
		result = prime * result + offset;
		return result;
	}
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		KmerSupportNode other = (KmerSupportNode) obj;
		if (evidence == null) {
			if (other.evidence != null)
				return false;
		} else if (!evidence.equals(other.evidence))
			return false;
		if (offset != other.offset)
			return false;
		return true;
	}
}