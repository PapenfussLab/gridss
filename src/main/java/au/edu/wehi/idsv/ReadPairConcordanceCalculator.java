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
	public abstract boolean isConcordant(SAMRecord local);
}
