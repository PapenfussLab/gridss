package au.edu.wehi.idsv.debruijn;

import gnu.trove.iterator.TLongIntIterator;
import gnu.trove.map.TLongIntMap;
import gnu.trove.map.hash.TLongIntHashMap;

import java.util.AbstractList;
import java.util.ArrayList;
import java.util.List;

import au.edu.wehi.idsv.BreakendSummary;
import au.edu.wehi.idsv.LinearGenomicCoordinate;
import au.edu.wehi.idsv.Models;
import au.edu.wehi.idsv.util.ArrayHelper;

import com.google.common.collect.Lists;
import com.google.common.collect.Ordering;
import com.google.common.primitives.Ints;

public class DeBruijnNodeBase {
	private long kmer;
	/**
	 * Cached total support for kmer
	 */
	private int cacheWeight;
	/**
	 * Cached total support for kmer
	 */
	private int cacheReferenceWeight;
	/**
	 * Contributing evidence
	 */
	private List<Support> support = new ArrayList<Support>(4);
	private static class Support {
		public Support(int weight, long expectedLinearPosition, String evidenceID, boolean isReference, int category, BreakendSummary breakend) {
			this.weight = weight;
			this.expectedLinearPosition = expectedLinearPosition;
			this.evidenceID = evidenceID;
			this.isReference = isReference;
			this.category = category;
			this.breakend = breakend;
		}
		public final long expectedLinearPosition;
		public final int weight;
		public final int category;
		public final boolean isReference;
		public final String evidenceID;
		public final BreakendSummary breakend;
	}
	public DeBruijnNodeBase(VariantEvidence evidence, int readKmerOffset, ReadKmer kmer) {
		this(kmer.kmer,
			kmer.weight,
			evidence.getExpectedLinearPosition(readKmerOffset),
			evidence.getEvidenceID(),
			evidence.isReferenceKmer(readKmerOffset),
			evidence.getCategory(),
			evidence.getBreakend());
	}
	/**
	 * Creates a node from the given read with the given level of support
	 * @param weight support weight
	 * @param read source read
	 */
	public DeBruijnNodeBase(long kmer, int weight, long expectedLinearPosition, String evidenceID, boolean isReference, int category, BreakendSummary breakend) {
		if (weight <= 0) throw new IllegalArgumentException("weight must be positive");
		this.support.add(new Support(weight, expectedLinearPosition, evidenceID, isReference, category, breakend));
		this.kmer = kmer;
		this.cacheWeight = weight;
		if (isReference) {
			this.cacheReferenceWeight = weight;
		}
	}
	/**
	 * Merges the given node into this one
	 * @param node
	 */
	public void add(DeBruijnNodeBase node) {
		this.support.addAll(node.support);
		this.cacheWeight += node.cacheWeight;
		this.cacheReferenceWeight += node.cacheReferenceWeight;
	}
	/**
	 * Reduces the weighting of this node due to removal of a supporting read
	 */
	public void remove(DeBruijnNodeBase node) {
		this.support.removeAll(node.support);
		this.cacheWeight -= node.cacheWeight;
		this.cacheReferenceWeight -= node.cacheReferenceWeight;
	}
	public long getKmer() {
		return kmer;
	}
	/**
	 * returns the weight of this node
	 * @return weight of this node
	 */
	public int getWeight() {
		return cacheWeight;
	}
	/**
	 * Indicates whether the given kmer provides any reference support
	 * @return true if kmer supports reference allele, false otherwise
	 */
	public boolean isReference() {
		return cacheReferenceWeight > 0;
	}
	/**
	 * Reads supporting this kmer. Reads containing this kmer multiple times will have multiple entries
	 * @return supporting reads
	 */
	public List<String> getSupportingEvidenceList() {
		return new SupportEvidenceIDList();
	}
	public static BreakendSummary getExpectedBreakend(LinearGenomicCoordinate lgc, Iterable<? extends DeBruijnNodeBase> path) {
		List<BreakendSummary> bs = Lists.newArrayList();
		List<Long> weights = Lists.newArrayList();
		for (DeBruijnNodeBase n : path) {
			bs.addAll(n.new SupportBreakendList());
			weights.addAll(n.new SupportWeightAsLongList());
		}
		return Models.calculateBreakend(lgc, bs, weights);
	}
	public long getExpectedPosition() {
		TLongIntMap map = new TLongIntHashMap();
		for (Support s : support) {
			map.adjustOrPutValue(s.expectedLinearPosition, s.weight, s.weight);
		}
		return getKeyWithMaxValue(map);
	}
	public long getExpectedReferencePosition() {
		if (!isReference()) throw new IllegalStateException("No reference kmer support");
		TLongIntMap map = new TLongIntHashMap();
		for (Support s : support) {
			if (s.isReference) {
				map.adjustOrPutValue(s.expectedLinearPosition, s.weight, s.weight);
			}
		}
		return getKeyWithMaxValue(map);
	}
	public static <T extends DeBruijnNodeBase> long getExpectedReferencePosition(Iterable<T> nodes) {
		TLongIntMap map = new TLongIntHashMap();
		for (DeBruijnNodeBase n : nodes) {
			for (Support s : n.support) {
				if (s.isReference) {
					map.adjustOrPutValue(s.expectedLinearPosition, s.weight, s.weight);
				}
			}
		}
		if (map.isEmpty()) throw new IllegalStateException("No reference support");
		return getKeyWithMaxValue(map);
	}
	private static long getKeyWithMaxValue(TLongIntMap map) {
		int maxValue = 0;
		long maxPos = 0;
		for (TLongIntIterator it = map.iterator(); it.hasNext(); ) {
		    it.advance();
		    if (it.value() > maxValue) {
		    	maxValue = it.value();
		    	maxPos = it.key();
		    }
		}
		return maxPos;
		
	}
	public int[] getCountByCategory() {
		int[] array = null;
		for (Support s : support) {
			array = ArrayHelper.add(array, s.category, 1);
		}
		return array;
	}
	public static Ordering<? extends DeBruijnNodeBase> ByWeight = new Ordering<DeBruijnNodeBase>() {
		public int compare(DeBruijnNodeBase o1, DeBruijnNodeBase o2) {
			  return Ints.compare(o1.getWeight(), o2.getWeight());
		  }
	};
	@Override
	public String toString() {
		String str = String.format("%s w=%d, #=%d p=%d",
				isReference() ? "R" : " ",
				cacheWeight,
				getSupportingEvidenceList().size(),
				getExpectedPosition());
		if (isReference() && getExpectedPosition() != getExpectedReferencePosition()) {
			str += String.format(" r=%d", getExpectedReferencePosition());
		}
		return str;
	}
	private class SupportEvidenceIDList extends AbstractList<String> {
		@Override
		public String get(int index) { return support.get(index).evidenceID; }
		@Override
		public int size() { return support.size(); }
	}
	private class SupportWeightAsLongList extends AbstractList<Long> {
		@Override
		public Long get(int index) { return (long)support.get(index).weight; }
		@Override
		public int size() { return support.size(); }
	}
	private class SupportBreakendList extends AbstractList<BreakendSummary> {
		@Override
		public BreakendSummary get(int index) { return support.get(index).breakend; }
		@Override
		public int size() { return support.size(); }
	}
}
