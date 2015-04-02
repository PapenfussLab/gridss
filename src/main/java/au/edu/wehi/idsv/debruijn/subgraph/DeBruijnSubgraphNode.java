package au.edu.wehi.idsv.debruijn.subgraph;

import au.edu.wehi.idsv.debruijn.DeBruijnNodeBase;
import au.edu.wehi.idsv.debruijn.ReadKmer;
import au.edu.wehi.idsv.debruijn.VariantEvidence;

public class DeBruijnSubgraphNode extends DeBruijnNodeBase {
	private long maxPosition;
	private SubgraphSummary subgraph;
	public DeBruijnSubgraphNode(VariantEvidence evidence, int readKmerOffset, ReadKmer kmer) {
		super(evidence, readKmerOffset, kmer);
		this.maxPosition = evidence.getExpectedLinearPosition(readKmerOffset);
	}
	/**
	 * Highest expected reference-allele genomic position
	 * @return
	 */
	public long getMaxLinearPosition() {
		return maxPosition;
	}
	@Override
	public void add(DeBruijnNodeBase node) {
		assert(node instanceof DeBruijnSubgraphNode);
		this.maxPosition = Math.max(maxPosition, ((DeBruijnSubgraphNode)node).maxPosition);
		super.add(node);
	}
	public boolean remove(DeBruijnSubgraphNode node) {
		throw new UnsupportedOperationException("Unable to remove read support from individual kmers: unable to calculate maxPosition.");
	}
	public SubgraphSummary getSubgraph() {
		if (subgraph == null) return null;
		return subgraph.getRoot();
	}
	public void setSubgraph(SubgraphSummary subgraph) {
		this.subgraph = subgraph;
	}
	public String toString() {
		return String.format("%s g=%d max=%d",
				super.toString(),
				subgraph.getAnyKmer(),
				maxPosition);
	}
}
