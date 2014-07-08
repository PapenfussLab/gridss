package au.edu.wehi.idsv.debruijn;

import htsjdk.samtools.SAMRecord;

import java.util.Set;

import com.google.common.collect.Sets;

public class DeBruijnNodeBase {
	private int nodeWeight = 0;
	// technically this is incorrect as we should have a MultiSet since
	// a record can contain a kmer multiple times but since we only add/remove
	// entire reads at a time, it works out ok for the current implementation
	private Set<SAMRecord> supportSet = Sets.newHashSet();
	public DeBruijnNodeBase(VariantEvidence evidence, int readKmerOffset, ReadKmer kmer) {
		if (kmer.weight <= 0) throw new IllegalArgumentException("weight must be positive");
		this.nodeWeight = kmer.weight;
		supportSet.add(evidence.getSAMRecord());
	}
	/**
	 * Merges the given node into this one
	 * @param node
	 */
	public void add(DeBruijnNodeBase node) {
		this.nodeWeight += node.nodeWeight;
		this.supportSet.addAll(node.supportSet);
	}
	/**
	 * Reduces the weighting of this node due to removal of a supporting read
	 * @return true if this node is now weightless
	 */
	public void remove(DeBruijnNodeBase node) {
		this.nodeWeight -= node.nodeWeight;
		this.supportSet.remove(node.supportSet);
	}
	/**
	 * returns the weight of this node
	 * @return weight of this node
	 */
	public int getWeight() {
		return nodeWeight;
	}
	/**
	 * Reads supporting this kmer
	 * @return supporting reads
	 */
	public Set<SAMRecord> getSupportingReads() {
		return supportSet;
	}
	@Override
	public String toString() {
		return String.format("w=%d, #reads=%d", nodeWeight, supportSet.size());
	}
}
