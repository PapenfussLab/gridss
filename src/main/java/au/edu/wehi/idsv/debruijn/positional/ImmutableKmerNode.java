package au.edu.wehi.idsv.debruijn.positional;

public class ImmutableKmerNode implements KmerNode {
	private final long kmer;
	private final int start;
	private final int end;
	private final int weight;
	private boolean isReference;
	public ImmutableKmerNode(
			long kmer,
			int start,
			int end,
			int weight,
			boolean isReference) {
		this.kmer = kmer;
		this.start = start;
		this.end = end;
		this.weight = weight;
		this.isReference = isReference;
	}
	public ImmutableKmerNode(KmerNode node) {
		this(node.kmer(), node.startPosition(), node.endPosition(), node.weight(), node.isReference());
	}
	public long kmer() { return kmer; }
	public int startPosition() { return start; }
	public int endPosition() { return end; }
	public int weight() { return weight; }
	public boolean isReference() { return isReference; }
}
