package au.edu.wehi.idsv.debruijn;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.BreakendSummary;
import au.edu.wehi.idsv.LinearGenomicCoordinate;
import au.edu.wehi.idsv.Models;
import au.edu.wehi.idsv.util.ArrayHelper;

import com.google.common.collect.Lists;
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
	private List<String> supportList = new ArrayList<String>(4);
	private long tPositionWeightSum = 0;
	private long fPositionWeightSum = 0;
	private long bPositionWeightSum = 0;
	private int tWeightSum = 0;
	private int fWeightSum = 0;
	private int bWeightSum = 0;
	private int[] categoryCount;
	private List<BreakendSummary> breakends = new ArrayList<BreakendSummary>(4);
	private List<Long> weights = new ArrayList<Long>(4);
	public DeBruijnNodeBase(VariantEvidence evidence, int readKmerOffset, ReadKmer kmer) {
		this(kmer.weight,
			evidence.getExpectedLinearPosition(readKmerOffset),
			evidence.getEvidenceID(),
			evidence.isReferenceKmer(readKmerOffset),
			evidence.isDirectlyAnchoredToReference() ? evidence.getBreakend().direction : null,
			evidence.getCategory(),
			evidence.getBreakend());
	}
	/**
	 * Creates a node from the given read with the given level of support
	 * @param weight support weight
	 * @param read source read
	 */
	public DeBruijnNodeBase(int weight, long expectedLinearPosition, String evidenceID, boolean isReference, BreakendDirection anchorDirection, int category, BreakendSummary unanchoredBreakend) {
		if (weight <= 0) throw new IllegalArgumentException("weight must be positive");
		assert(unanchoredBreakend != null);
		this.nodeWeight = weight;
		supportList.add(evidenceID);
		if (isReference) referenceSupport++;
		tPositionWeightSum += expectedLinearPosition * weight;
		tWeightSum += weight;
		if (anchorDirection == BreakendDirection.Forward) {
			fPositionWeightSum += expectedLinearPosition * weight;
			fWeightSum += weight;
		} else if (anchorDirection == BreakendDirection.Backward) {
			bPositionWeightSum += expectedLinearPosition * weight;
			bWeightSum += weight;
		}
		categoryCount = new int[category + 1];
		categoryCount[category] = 1;
		breakends.add(unanchoredBreakend);
		weights.add((long)weight);
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
		this.categoryCount = ArrayHelper.add(this.categoryCount, node.categoryCount);
		this.breakends.addAll(node.breakends);
		this.weights.addAll(node.weights);
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
		this.categoryCount = ArrayHelper.subtract(this.categoryCount, node.categoryCount);
		this.breakends.removeAll(node.breakends);
		this.weights.removeAll(node.weights);
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
	public List<String> getSupportingEvidenceList() {
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
	public static BreakendSummary getExpectedBreakend(LinearGenomicCoordinate lgc, Collection<? extends DeBruijnNodeBase> path) {
		List<BreakendSummary> bs = Lists.newArrayList();
		List<Long> weights = Lists.newArrayList();
		for (DeBruijnNodeBase n : path) {
			bs.addAll(n.breakends);
			weights.addAll(n.weights);
		}
		return Models.calculateBreakend(lgc, bs, weights);
	}
	public double getExpectedPosition() {
		return tPositionWeightSum / (double)tWeightSum;
	}
	public int[] getCountByCategory() {
		return categoryCount;
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
