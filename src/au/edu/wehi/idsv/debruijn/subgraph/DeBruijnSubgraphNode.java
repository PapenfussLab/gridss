package au.edu.wehi.idsv.debruijn.subgraph;

import au.edu.wehi.idsv.debruijn.DeBruijnNodeBase;
import au.edu.wehi.idsv.debruijn.ReadKmer;
import au.edu.wehi.idsv.debruijn.VariantEvidence;

import com.google.common.collect.Ordering;
import com.google.common.primitives.Ints;

public class DeBruijnSubgraphNode extends DeBruijnNodeBase {
	/**
	 * Genomic coordinates of positions this kmer is aligned to the reference
	 * 
	 * Note: reads with less than k bases mapped to the reference will be considered unanchored 
	 */
	private int minReferencePosition = Integer.MAX_VALUE;
	private int maxReferencePosition = Integer.MIN_VALUE;
	/**
	 * Genomic coordinates of closest anchored mate of reads starting with this kmer 
	 */
	private int minMatePosition = Integer.MAX_VALUE;
	private int maxMatePosition = Integer.MIN_VALUE;
	 
	private SubgraphSummary subgraph;
	public DeBruijnSubgraphNode(VariantEvidence evidence, int readKmerOffset, ReadKmer kmer) {
		super(evidence, readKmerOffset, kmer);
		if (evidence.isReferenceKmer(readKmerOffset)) {
			maxReferencePosition = minReferencePosition = evidence.getInferredReferencePosition(readKmerOffset);
		} else if (!evidence.isDirectlyAnchoredToReference()) {
			minMatePosition = maxMatePosition = evidence.getMateAnchorPosition();
		}
	}
	public void add(DeBruijnSubgraphNode node) {
		super.add(node);
		minReferencePosition = Math.min(minReferencePosition, node.minReferencePosition);
		maxReferencePosition = Math.max(maxReferencePosition, node.maxReferencePosition);
		minMatePosition = Math.min(minMatePosition, node.minMatePosition);
		maxMatePosition = Math.max(maxMatePosition, node.maxMatePosition);
	}
	public boolean remove(DeBruijnSubgraphNode node) {
		throw new UnsupportedOperationException("Unable to remove read support from individual kmers: only most closest anchors are tracked.");
	}
	public boolean isReference() {
		return minReferencePosition != Integer.MAX_VALUE;
	}
	public boolean isMateAnchored() {
		return minMatePosition != Integer.MAX_VALUE;
	}
	public SubgraphSummary getSubgraph() {
		if (subgraph == null) return null;
		return subgraph.getRoot();
	}
	public void setSubgraph(SubgraphSummary subgraph) {
		this.subgraph = subgraph;
	}
	public Integer getMinReferencePosition() {
		return isReference() ? minReferencePosition : null;
	}
	public Integer getMaxReferencePosition() {
		return isReference() ? maxReferencePosition : null;
	}
	public Integer getMinMatePosition() {
		return isMateAnchored() ? minMatePosition : null;
	}
	public Integer getMaxMatePosition() {
		return isMateAnchored() ? maxMatePosition : null; 
	}
	public static Ordering<DeBruijnSubgraphNode> ByWeight = new Ordering<DeBruijnSubgraphNode>() {
		public int compare(DeBruijnSubgraphNode o1, DeBruijnSubgraphNode o2) {
			  return Ints.compare(o1.getWeight(), o2.getWeight());
		  }
	};
	public String toString() {
		return String.format("%s%s g=%d,%s",
				maxReferencePosition == Integer.MIN_VALUE ? " " : "R",
				maxMatePosition == Integer.MIN_VALUE ? " " : "M",
				subgraph.getAnyKmer(),
				super.toString());
	}
}
