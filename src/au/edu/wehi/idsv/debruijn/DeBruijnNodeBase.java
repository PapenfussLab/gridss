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
	/**
	 * Genomic coordinates of positions this kmer is aligned to the reference
	 * 
	 * Note: reads with less than k bases mapped to the reference will be considered unanchored 
	 */
	private SortedMultiset<Integer> referenceCoordinateSet;
	/**
	 * Genomic coordinates of unanchored reads containing (TODO: or just starting with?) this kmer 
	 */
	private SortedMultiset<Integer> unanchoredReferenceCoordinateSet;
	public void add(DeBruijnEvidence evidence, int readKmerOffset, ReadKmer kmer) {
		if (kmer.weight <= 0) throw new IllegalArgumentException("weight must be positive");
		nodeWeight += kmer.weight;
		supportSet.add(evidence.getSAMRecord());
		if (evidence.isReferenceAnchor(readKmerOffset)) referenceCoordinateSet.add(evidence.getReferenceKmerAnchorPosition());
		if (evidence.isMateAnchor(readKmerOffset)) unanchoredReferenceCoordinateSet.add(evidence.getMateAnchorPosition());
	}
	/**
	 * Reduces the weighting of this node due to removal of a supporting read
	 * @return true if this node is now weightless
	 */
	public boolean remove(DeBruijnEvidence evidence, int readKmerOffset, ReadKmer kmer) {
		if (kmer.weight <= 0) throw new IllegalArgumentException("weight must be positive");
		nodeWeight -= kmer.weight;
		supportSet.remove(evidence.getSAMRecord());
		if (evidence.isReferenceAnchor(readKmerOffset)) referenceCoordinateSet.remove(evidence.getReferenceKmerAnchorPosition());
		if (evidence.isMateAnchor(readKmerOffset)) unanchoredReferenceCoordinateSet.remove(evidence.getMateAnchorPosition());
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
