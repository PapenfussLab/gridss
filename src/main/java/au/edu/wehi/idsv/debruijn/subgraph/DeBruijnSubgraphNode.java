package au.edu.wehi.idsv.debruijn.subgraph;

import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TIntIntHashMap;
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
	TIntList referencePosition = new TIntArrayList(4);
	TIntList referencePositionWeight = new TIntArrayList(4);
	private int minReferencePosition = Integer.MAX_VALUE;
	private int maxReferencePosition = Integer.MIN_VALUE;
	/**
	 * Genomic coordinates of closest anchored mate of reads starting with this kmer 
	 */
	//private TIntList matePosition = new TIntArrayList();
	//private TIntList matePositionWeight = new TIntArrayList();
	private int minMatePosition = Integer.MAX_VALUE;
	private int maxMatePosition = Integer.MIN_VALUE;
	 
	private SubgraphSummary subgraph;
	public DeBruijnSubgraphNode(VariantEvidence evidence, int readKmerOffset, ReadKmer kmer) {
		super(evidence, readKmerOffset, kmer);
		if (evidence.isReferenceKmer(readKmerOffset)) {
			addRefPosition(evidence.getInferredReferencePosition(readKmerOffset), kmer.weight);
		} else if (!evidence.isDirectlyAnchoredToReference()) {
			addMatePosition(evidence.getMateAnchorPosition(), kmer.weight);
		}
	}
	private void addMatePosition(int position, int weight) {
		minMatePosition = Math.min(minMatePosition, position);
		maxMatePosition = Math.max(maxMatePosition, position);
		//matePosition.add(position);
		//matePositionWeight.add(weight);
	}
	private void addRefPosition(int position, int weight) {
		minReferencePosition = Math.min(minReferencePosition, position);
		maxReferencePosition = Math.max(maxReferencePosition, position);
		referencePosition.add(position);
		referencePositionWeight.add(weight);
	}
	private static Integer getBestPosition(TIntList positions, TIntList weights) {
		TIntIntHashMap lookup = new TIntIntHashMap();
		int bestWeight = 0;
		int bestPos = -1;
		for (int i = 0; i < positions.size(); i++) {
			int added = lookup.adjustOrPutValue(positions.get(i), weights.get(i), weights.get(i));
			if (added > bestWeight) {
				bestWeight = added;
				bestPos = positions.get(i);
			}
		}
		if (bestWeight > 0) return bestPos;
		return null;
	}
	@Override
	public void add(DeBruijnNodeBase node) {
		assert(node instanceof DeBruijnSubgraphNode);
		super.add(node);
		DeBruijnSubgraphNode s = (DeBruijnSubgraphNode)node;
		this.referencePosition.addAll(s.referencePosition);
		this.referencePositionWeight.addAll(s.referencePositionWeight);
		//this.matePosition.addAll(s.matePosition);
		//this.matePositionWeight.addAll(s.matePositionWeight);
		this.minReferencePosition = Math.min(minReferencePosition, s.minReferencePosition);
		this.maxReferencePosition = Math.max(maxReferencePosition, s.maxReferencePosition);
		this.minMatePosition = Math.min(minMatePosition, s.minMatePosition);
		this.maxMatePosition = Math.max(maxMatePosition, s.maxMatePosition);
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
	public Integer getBestReferencePosition() {
		return getBestPosition(referencePosition, referencePositionWeight);
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
				isReference() ? "R" : " ",
				isMateAnchored() ? "M" : " ",
				subgraph.getAnyKmer(),
				super.toString());
	}
}
