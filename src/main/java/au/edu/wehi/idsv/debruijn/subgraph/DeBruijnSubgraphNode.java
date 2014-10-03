package au.edu.wehi.idsv.debruijn.subgraph;

import java.util.NavigableMap;

import au.edu.wehi.idsv.debruijn.DeBruijnNodeBase;
import au.edu.wehi.idsv.debruijn.ReadKmer;
import au.edu.wehi.idsv.debruijn.VariantEvidence;

import com.google.common.collect.Maps;
import com.google.common.collect.Ordering;
import com.google.common.primitives.Ints;

public class DeBruijnSubgraphNode extends DeBruijnNodeBase {
	/**
	 * Genomic coordinates of positions this kmer is aligned to the reference
	 * 
	 * Note: reads with less than k bases mapped to the reference will be considered unanchored 
	 */
	private NavigableMap<Integer, Integer> referencePosition = Maps.newTreeMap();
	/**
	 * Genomic coordinates of closest anchored mate of reads starting with this kmer 
	 */
	private NavigableMap<Integer, Integer> matePosition = Maps.newTreeMap();
	 
	private SubgraphSummary subgraph;
	public DeBruijnSubgraphNode(VariantEvidence evidence, int readKmerOffset, ReadKmer kmer) {
		super(evidence, readKmerOffset, kmer);
		if (evidence.isReferenceKmer(readKmerOffset)) {
			addPosition(referencePosition, evidence.getInferredReferencePosition(readKmerOffset), kmer.weight);
		} else if (!evidence.isDirectlyAnchoredToReference()) {
			addPosition(matePosition, evidence.getMateAnchorPosition(), kmer.weight);
		}
	}
	private static void addPosition(NavigableMap<Integer, Integer> positionMap, int position, int weight) {
		Integer existing = positionMap.get(position);
		if (existing == null) existing = 0;
		positionMap.put(position, existing + weight);
	}
	private static void addPositions(NavigableMap<Integer, Integer> positionMap, NavigableMap<Integer, Integer> toAdd) {
		for (int key : toAdd.keySet()) {
			addPosition(positionMap, key, toAdd.get(key));
		}
	}
	private static Integer getBestPosition(NavigableMap<Integer, Integer> positionMap) {
		int bestWeight = 0;
		int bestPos = -1;
		for (int pos : positionMap.keySet()) {
			int weight = positionMap.get(pos); 
			if (weight >= bestWeight) {
				bestPos = pos;
				bestWeight = weight;
			}
		}
		if (bestWeight > 0) return bestPos;
		return null;
	}
	@Override
	public void add(DeBruijnNodeBase node) {
		assert(node instanceof DeBruijnSubgraphNode);
		super.add(node);
		addPositions(referencePosition, ((DeBruijnSubgraphNode)node).referencePosition);
		addPositions(matePosition, ((DeBruijnSubgraphNode)node).matePosition);
	}
	public boolean remove(DeBruijnSubgraphNode node) {
		throw new UnsupportedOperationException("Unable to remove read support from individual kmers: only most closest anchors are tracked.");
	}
	public boolean isReference() {
		return !referencePosition.isEmpty();
	}
	public boolean isMateAnchored() {
		return !matePosition.isEmpty();
	}
	public SubgraphSummary getSubgraph() {
		if (subgraph == null) return null;
		return subgraph.getRoot();
	}
	public void setSubgraph(SubgraphSummary subgraph) {
		this.subgraph = subgraph;
	}
	public Integer getMinReferencePosition() {
		return isReference() ? referencePosition.firstKey() : null;
	}
	public Integer getMaxReferencePosition() {
		return isReference() ? referencePosition.lastKey() : null;
	}
	public Integer getBestReferencePosition() {
		return getBestPosition(referencePosition);
	}
	public Integer getMinMatePosition() {
		return isMateAnchored() ? matePosition.firstKey() : null;
	}
	public Integer getMaxMatePosition() {
		return isMateAnchored() ? matePosition.lastKey() : null; 
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
