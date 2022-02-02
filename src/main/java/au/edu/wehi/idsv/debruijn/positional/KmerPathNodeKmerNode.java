package au.edu.wehi.idsv.debruijn.positional;

/**
 * KmerNode representing an individual kmer of a KmerPathNode
 * @author Daniel Cameron
 *
 */
public class KmerPathNodeKmerNode implements KmerNode {
	private final KmerPathNode node;
	/**
	 * Offset indicates which kmer in the KmerPathNode this kmer corresponds to
	 * 
	 * A negative value indicates a collapsed kmer
	 */
	private final int offset;
	public KmerPathNodeKmerNode(KmerPathNode node, int offset) {
		this.node = node;
		this.offset = offset;
		assert(offset < node.length());
	}
	public KmerPathNode node() {
		return node;
	}
	public int offset() { return offset; }
	@Override
	public long lastKmer() {
		return node.kmer(offset);
	}
	@Override
	public int lastStart() {
		return node.startPosition(offset);
	}
	@Override
	public int lastEnd() {
		return node.endPosition(offset);
	}
	@Override
	public int weight() {
		return node.weight(offset);
	}
	@Override
	public boolean isReference() {
		return node.isReference();
	}
	@Override
	public String toString() {
		return String.format("{%d}%s", offset, node);
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		long kmer = node.firstKmer();
		int result = (int) (kmer ^ (kmer >>> 32));
		result = prime * result + node.firstStart();
		result = prime * result + offset;
		return result;
	}
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		KmerPathNodeKmerNode other = (KmerPathNodeKmerNode) obj;
		if (other == null)
			return false;
		return offset == other.offset && node == other.node; 
	}
}
