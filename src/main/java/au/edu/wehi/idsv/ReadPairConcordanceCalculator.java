package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;

/**
 * Determines whether the given read pair is considered concordant
 * @author Daniel Cameron
 *
 */
public abstract class ReadPairConcordanceCalculator {
	/**
	 * Fragment size of max concordant fragment
	 * @return
	 */
	public abstract int maxConcordantFragmentSize();
	public abstract int minConcordantFragmentSize();
	public boolean isConcordant(SAMRecord read) {
		return isConcordant(read, null);
	}
	public boolean isConcordant(SAMRecord read1, SAMRecord read2) {
		// (assumes FR)
		if (read2 == null) {
			return read1.getReadPairedFlag()
					&& !read1.getReadUnmappedFlag()
					&& !read1.getMateUnmappedFlag()
					&& read1.getReferenceIndex() == read1.getMateReferenceIndex()
					&& read1.getReadNegativeStrandFlag() != read1.getMateNegativeStrandFlag();
		} else {
			return read1.getReadPairedFlag()
					&& !read1.getReadUnmappedFlag()
					&& !read2.getReadUnmappedFlag()
					&& read1.getReferenceIndex() == read2.getReferenceIndex()
					&& read1.getReadNegativeStrandFlag() != read2.getReadNegativeStrandFlag(); 
		}
	}
}
