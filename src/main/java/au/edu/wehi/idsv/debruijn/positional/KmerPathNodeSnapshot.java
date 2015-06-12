package au.edu.wehi.idsv.debruijn.positional;


/**
 * Snapshot of a kmer path node used for tracking mutated nodes
 * 
 * @author cameron.d
 *
 */
class KmerPathNodeSnapshot extends ImmutableKmerNode {
	public final int version;
	public final KmerPathNode node;
	public KmerPathNodeSnapshot(KmerPathNode node) {
		super(node);
		this.node = node;
		this.version = node.version();
	}
	/**
	 * Node has remained unchanged since this snapshot was created
	 * @return
	 */
	public boolean isUnchanged() {
		return version == node.version();
	}
}