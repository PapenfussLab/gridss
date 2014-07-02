package au.edu.wehi.idsv.debruijn;

import htsjdk.samtools.SAMRecord;

import java.util.Set;

import com.google.common.collect.Sets;
import com.google.common.collect.SortedMultiset;

public class DeBruijnNodeBase {
	private long nodeWeight = 0;
	// technically this is incorrect as we should have a MultiSet since
	// a record can contain a kmer multiple times but since we only add/remove
	// entire reads at a time, it works out ok for the current implementation
	private Set<SAMRecord> supportSet = Sets.newHashSet();
	public void add(DeBruijnEvidence evidence, int readKmerOffset, ReadKmer kmer) {
		if (kmer.weight <= 0) throw new IllegalArgumentException("weight must be positive");
		nodeWeight += kmer.weight;
		supportSet.add(evidence.getSAMRecord());
	}
	/**
	 * Reduces the weighting of this node due to removal of a supporting read
	 * @return true if this node is now weightless
	 */
	public boolean remove(DeBruijnEvidence evidence, int readKmerOffset, ReadKmer kmer) {
		if (kmer.weight <= 0) throw new IllegalArgumentException("weight must be positive");
		nodeWeight -= kmer.weight;
		supportSet.remove(evidence.getSAMRecord());
		return nodeWeight <= 0;
	}
	/**
	 * returns the weight of this node
	 * @return weight of this node
	 */
	public long getWeight() {
		return nodeWeight;
	}
	/**
	 * Reads supporting this kmer
	 * @return supporting reads
	 */
	public Set<SAMRecord> getSupportingReads() {
		return supportSet;
	}
}
