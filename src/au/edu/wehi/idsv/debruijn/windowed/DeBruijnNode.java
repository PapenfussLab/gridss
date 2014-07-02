package au.edu.wehi.idsv.debruijn.windowed;

import au.edu.wehi.idsv.debruijn.DeBruijnEvidence;
import au.edu.wehi.idsv.debruijn.DeBruijnNodeBase;
import au.edu.wehi.idsv.debruijn.ReadKmer;

import com.google.common.collect.SortedMultiset;

public class DeBruijnNode extends DeBruijnNodeBase {
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
	@Override
	public void add(DeBruijnEvidence evidence, int readKmerOffset, ReadKmer kmer) {
		super.add(evidence, readKmerOffset, kmer);
		if (evidence.isReferenceAnchor(readKmerOffset)) referenceCoordinateSet.add(evidence.getReferenceKmerAnchorPosition());
		if (evidence.isMateAnchor(readKmerOffset)) unanchoredReferenceCoordinateSet.add(evidence.getMateAnchorPosition());
	}
	/**
	 * Reduces the weighting of this node due to removal of a supporting read
	 * @return true if this node is now weightless
	 */
	@Override
	public boolean remove(DeBruijnEvidence evidence, int readKmerOffset, ReadKmer kmer) {
		if (evidence.isReferenceAnchor(readKmerOffset)) referenceCoordinateSet.remove(evidence.getReferenceKmerAnchorPosition());
		if (evidence.isMateAnchor(readKmerOffset)) unanchoredReferenceCoordinateSet.remove(evidence.getMateAnchorPosition());
		return super.remove(evidence, readKmerOffset, kmer);
	}
}
