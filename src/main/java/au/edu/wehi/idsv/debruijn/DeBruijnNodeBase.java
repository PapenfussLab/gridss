package au.edu.wehi.idsv.debruijn;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.DirectedEvidence;

import com.google.common.collect.Ordering;
import com.google.common.primitives.Ints;

public class DeBruijnNodeBase {
	/**
	 * Total support for kmer
	 */
	private int nodeWeight = 0;
	/**
	 * Reference kmer support
	 */
	private int referenceSupport = 0;
	/**
	 * Contributing evidence
	 */
	private List<DirectedEvidence> supportList = new ArrayList<DirectedEvidence>(4);
	private long tPositionWeightSum = 0;
	private long fPositionWeightSum = 0;
	private long bPositionWeightSum = 0;
	private int tWeightSum = 0;
	private int fWeightSum = 0;
	private int bWeightSum = 0;
	public DeBruijnNodeBase(VariantEvidence evidence, int readKmerOffset, ReadKmer kmer) {
		this(kmer.weight, evidence.getExpectedLinearPosition(readKmerOffset), evidence.getDirectedEvidence(), evidence.isReferenceKmer(readKmerOffset), evidence.isDirectlyAnchoredToReference());
	}
	/**
	 * Creates a node from the given read with the given level of support
	 * @param weight support weight
	 * @param read source read
	 */
	public DeBruijnNodeBase(int weight, long expectedLinearPosition, DirectedEvidence evidence, boolean isReference, boolean isAnchored) {
		if (weight <= 0) throw new IllegalArgumentException("weight must be positive");
		this.nodeWeight = weight;
		supportList.add(evidence);
		if (isReference) referenceSupport++;
		tPositionWeightSum += expectedLinearPosition * weight;
		tWeightSum += weight;
		if (isAnchored) {
			if (evidence.getBreakendSummary().direction == BreakendDirection.Forward) {
				fPositionWeightSum += expectedLinearPosition * weight;
				fWeightSum += weight;
			} else {
				bPositionWeightSum += expectedLinearPosition * weight;
				bWeightSum += weight;
			}
		}
	}
	/**
	 * Merges the given node into this one
	 * @param node
	 */
	public void add(DeBruijnNodeBase node) {
		this.nodeWeight += node.nodeWeight;
		this.referenceSupport += node.referenceSupport;
		this.supportList.addAll(node.supportList);
		this.tPositionWeightSum += node.tPositionWeightSum;
		this.fPositionWeightSum += node.fPositionWeightSum;
		this.bPositionWeightSum += node.bPositionWeightSum;
		this.tWeightSum += node.tWeightSum;
		this.fWeightSum += node.fWeightSum;
		this.bWeightSum += node.bWeightSum;
	}
	/**
	 * Reduces the weighting of this node due to removal of a supporting read
	 */
	public void remove(DeBruijnNodeBase node) {
		this.nodeWeight -= node.nodeWeight;
		this.referenceSupport -= node.referenceSupport;
		this.supportList.remove(node.supportList);
		this.tPositionWeightSum -= node.tPositionWeightSum;
		this.fPositionWeightSum -= node.fPositionWeightSum;
		this.bPositionWeightSum -= node.bPositionWeightSum;
		this.tWeightSum -= node.tWeightSum;
		this.fWeightSum -= node.fWeightSum;
		this.bWeightSum -= node.bWeightSum;
	}
	/**
	 * returns the weight of this node
	 * @return weight of this node
	 */
	public int getWeight() {
		return nodeWeight;
	}
	/**
	 * Indicates whether the given kmer provides any reference support
	 * @return true if kmer supports reference allele, false otherwise
	 */
	public boolean isReference() {
		return referenceSupport > 0;
	}
	/**
	 * Reads supporting this kmer. Reads containing this kmer multiple times will have multiple entries
	 * @return supporting reads
	 */
	public List<DirectedEvidence> getSupportingEvidenceList() {
		return supportList;
	}
	/**
	 * Gets the inferred reference anchor position of the base immediately adjacent to the given breakend sequence
	 * @param direction breakend direction.
	 * Forward returns the inferred position of the kmer immediately before the start of the path
	 * Backwards returns the inferred position of the kmer immediately after the end of the path
	 * @param path kmer path
	 * @return Inferred reference anchor position 
	 */
	public static double getExpectedPositionForDirectAnchor(BreakendDirection direction, Collection<? extends DeBruijnNodeBase> path) {
		long posWeightSum = 0;
		long weightSum = 0;
		//long oppositePosWeightSum = 0;
		//long oppositeWeightSum = 0;
		int offset;
		if (direction == BreakendDirection.Forward) {
			offset = 1;
			for (DeBruijnNodeBase n : path) {
				weightSum += n.fWeightSum;
				posWeightSum += n.fPositionWeightSum - n.fWeightSum * offset;
				//oppositeWeightSum += n.bWeightSum;
				//oppositePosWeightSum += n.bPositionWeightSum - n.bWeightSum * offset;
				offset++;
			}
		} else {
			offset = -path.size();
			for (DeBruijnNodeBase n : path) {
				weightSum += n.bWeightSum;
				posWeightSum += n.bPositionWeightSum - n.bWeightSum * offset;
				//oppositeWeightSum += n.fWeightSum;
				//oppositePosWeightSum += n.fPositionWeightSum - n.fWeightSum * offset;
				offset++;
			}
		}
		double result = posWeightSum / (double)weightSum;
		return result;
	}
	public double getExpectedAnchorPosition() {
		return (fPositionWeightSum + bPositionWeightSum) / (double)(fWeightSum + bWeightSum); 
	}
	public double getExpectedPosition() {
		return tPositionWeightSum / (double)tWeightSum;
	}
	public static Ordering<? extends DeBruijnNodeBase> ByWeight = new Ordering<DeBruijnNodeBase>() {
		public int compare(DeBruijnNodeBase o1, DeBruijnNodeBase o2) {
			  return Ints.compare(o1.getWeight(), o2.getWeight());
		  }
	};
	@Override
	public String toString() {
		return String.format("%s w=%d, #=%d f=%.1f,b=%.1f",
				isReference() ? "R" : " ",
				nodeWeight,
				getSupportingEvidenceList().size(),
				fPositionWeightSum / (double)fWeightSum,
				bPositionWeightSum / (double)bWeightSum);
	}
}
