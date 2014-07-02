package au.edu.wehi.idsv.debruijn.windowed;

import com.google.common.collect.SortedMultiset;

import htsjdk.samtools.SAMRecord;

public class DeBruijnNode {
	//public boolean direction; // implicit in the parent graphs
	// store at node level? contig?
	/**
	 * Genomic coordinates of positions this kmer is aligned to the reference
	 * TODO: do we need offsets for these anchors as well? 
	 */
	public SortedMultiset<Integer> referenceCoordinate;
	/**
	 * Genomic coordinates of unanchored reads containing (TODO: or just starting with?) this kmer 
	 */
	public SortedMultiset<Integer> unanchoredReferenceCoordinate;
	/**
	 * Kmer weight
	 */
	public int weight;
	
}
