package au.edu.wehi.idsv.debruijn.subgraph;

import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.debruijn.DeBruijnNodeBase;
import au.edu.wehi.idsv.debruijn.ReadKmer;
import au.edu.wehi.idsv.debruijn.VariantEvidence;

import com.google.common.collect.Ordering;
import com.google.common.primitives.Ints;

public class DeBruijnSubgraphNode extends DeBruijnNodeBase {
	/**
	 * Genomic coordinate the positions this kmer would be aligned to
	 * 
	 */
	private TIntList fposition = new TIntArrayList(4);
	private TIntList fweight = new TIntArrayList(4);
	private TIntList bposition = new TIntArrayList(4);
	private TIntList bweight = new TIntArrayList(4);
	private int maxPosition;
	private int referenceSupport = 0;
	 
	private SubgraphSummary subgraph;
	public DeBruijnSubgraphNode(VariantEvidence evidence, int readKmerOffset, ReadKmer kmer) {
		super(evidence, readKmerOffset, kmer);
		int expectedPos = evidence.getExpectedReferencePosition(readKmerOffset);
		maxPosition = expectedPos;
		// TODO: split f and b evidence
		if (evidence.getDirection() == BreakendDirection.Forward) {
			fposition.add(expectedPos);
			fweight.add(kmer.weight);
		} else {
			bposition.add(expectedPos);
			bweight.add(kmer.weight);
		}
		if (evidence.isReferenceKmer(readKmerOffset)) {
			referenceSupport++;
		}
	}
	public int getMaxPosition() {
		return maxPosition;
	}
	@Override
	public void add(DeBruijnNodeBase node) {
		assert(node instanceof DeBruijnSubgraphNode);
		super.add(node);
		DeBruijnSubgraphNode s = (DeBruijnSubgraphNode)node;
		this.fposition.addAll(s.fposition);
		this.fweight.addAll(s.fweight);
		this.bposition.addAll(s.bposition);
		this.bweight.addAll(s.bweight);
		this.maxPosition = Math.max(maxPosition, s.maxPosition);
		this.referenceSupport += s.referenceSupport;
	}
	public boolean remove(DeBruijnSubgraphNode node) {
		throw new UnsupportedOperationException("Unable to remove read support from individual kmers: only most closest anchors are tracked.");
	}
	public boolean isReference() {
		return referenceSupport > 0;
	}
	public SubgraphSummary getSubgraph() {
		if (subgraph == null) return null;
		return subgraph.getRoot();
	}
	public void setSubgraph(SubgraphSummary subgraph) {
		this.subgraph = subgraph;
	}
	public static Ordering<DeBruijnSubgraphNode> ByWeight = new Ordering<DeBruijnSubgraphNode>() {
		public int compare(DeBruijnSubgraphNode o1, DeBruijnSubgraphNode o2) {
			  return Ints.compare(o1.getWeight(), o2.getWeight());
		  }
	};
	private static long weightedSum(TIntList value, TIntList weight, int offset) {
		long sum = 0;
		for (int i = 0; i < value.size(); i++) {
			long currentPos = value.get(i) + offset;
			long currentWeight = weight.get(i);
			sum += currentPos * currentWeight;
		}
		return sum;
	}
	public String toString() {
		return String.format("%s%s g=%d,f=%.1f,b=%.1f",
				isReference() ? "R" : " ",
				weightedSum(fposition, fweight, 0) / (double)fweight.sum(),
				weightedSum(bposition, bweight, 0) / (double)bweight.sum(),
				subgraph.getAnyKmer(),
				super.toString());
	}
}
