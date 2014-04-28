package au.edu.wehi.socrates.debruijn;

import java.util.Set;

import net.sf.samtools.SAMRecord;

import com.google.common.collect.Sets;

public class DeBruijnNode {
	private long nodeWeight = 0;
	// technically this is incorrect as we should have a MultiSet since
	// a record can contain a kmer multiple times but since we only add/remove
	// entire reads at a time, it works out ok for the current implementation
	private Set<SAMRecord> support = Sets.newHashSet();
	/**
	 * Adds weight to this node due to support of an additional read
	 */
	public void addRead(long weight, SAMRecord record) {
		if (weight <= 0) throw new IllegalArgumentException("weight must be positive");
		nodeWeight += weight;
		support.add(record);
	}
	/**
	 * Reduces the weighting of this node due to removal of a supporting read
	 * @return true if this node is now weightless
	 */
	public boolean removeRead(long weight, SAMRecord record) {
		if (weight <= 0) throw new IllegalArgumentException("weight must be positive");
		nodeWeight -= weight;
		support.remove(record);
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
		return support;
	}
}
