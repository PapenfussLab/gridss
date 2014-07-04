package au.edu.wehi.idsv.debruijn.subgraph;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.debruijn.DeBruijnEvidence;
import au.edu.wehi.idsv.debruijn.DeBruijnNodeBase;
import au.edu.wehi.idsv.debruijn.ReadKmer;

import com.google.common.collect.Ordering;
import com.google.common.primitives.Ints;

public class DeBruijnNode extends DeBruijnNodeBase {
	/**
	 * Genomic coordinates of positions this kmer is aligned to the reference
	 * 
	 * Note: reads with less than k bases mapped to the reference will be considered unanchored 
	 */
	private Integer referencePosition;
	/**
	 * Genomic coordinates of closest anchored mate of reads starting with this kmer 
	 */
	private Integer matePosition;
	private SubgraphSummary subgraph;
	@Override
	public void add(DeBruijnEvidence evidence, int readKmerOffset, ReadKmer kmer) {
		super.add(evidence, readKmerOffset, kmer);
		if (evidence.getDirection() == BreakendDirection.Forward) {
			// set to max anchor position
			if (evidence.isReferenceKmer(readKmerOffset)) {
				referencePosition = closest(evidence.getDirection(), referencePosition, evidence.getInferredReferencePosition(readKmerOffset));
			} else if (!evidence.isDirectlyAnchoredToReference()) {
				matePosition = closest(evidence.getDirection(), matePosition, evidence.getMateAnchorPosition());
			}
		}
	}
	private int closest(BreakendDirection direction, Integer current, int update) {
		if (current == null) return update;
		if (direction == BreakendDirection.Forward) {
			return Math.max(current, update);
		} else {
			return Math.min(current, update);
		}
	}
	/**
	 * Reduces the weighting of this node due to removal of a supporting read
	 * @return true if this node is now weightless
	 */
	@Override
	public boolean remove(DeBruijnEvidence evidence, int readKmerOffset, ReadKmer kmer) {
		throw new UnsupportedOperationException("Unable to remove read support from individual kmers: only most closest anchors are tracked.");
	}
	public boolean isReference() {
		return referencePosition != null;
	}
	public SubgraphSummary getSubgraph() {
		return subgraph.getRoot();
	}
	public void setSubgraph(SubgraphSummary subgraph) {
		this.subgraph = subgraph;
	}
	public Integer getReferencePosition() {
		return referencePosition;
	}
	public Integer getMatePosition() {
		return matePosition;
	}
	public static Ordering<DeBruijnNode> ByWeight = new Ordering<DeBruijnNode>() {
		public int compare(DeBruijnNode o1, DeBruijnNode o2) {
			  return Ints.compare(o1.getWeight(), o2.getWeight());
		  }
	};
}
