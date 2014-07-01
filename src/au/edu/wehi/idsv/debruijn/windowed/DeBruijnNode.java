package au.edu.wehi.idsv.debruijn.windowed;

import htsjdk.samtools.SAMRecord;

public class DeBruijnNode {
	//public boolean direction; // implicit in the parent graphs
	// store at node level? contig?
	/**
	 * Genomic coordinates of positions this kmer is aligned to the reference
	 * TODO: do we need offsets for these anchors as well? 
	 */
	public SortedMultiSet<Integer> referenceCoordinate;
	/**
	 * Genomic coordinates of unanchored reads containing (TODO: or just starting with?) this kmer 
	 */
	public SortedMultiSet<Integer> unanchoredReferenceCoordinate;
	/**
	 * Kmer weight
	 */
	public int weight;
	
}
