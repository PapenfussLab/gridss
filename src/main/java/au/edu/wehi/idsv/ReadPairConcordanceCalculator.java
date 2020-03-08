package au.edu.wehi.idsv;

import gridss.analysis.IdsvMetrics;
import gridss.analysis.InsertSizeDistribution;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Log;

/**
 * Determines whether the given read pair is considered concordant
 * @author Daniel Cameron
 *
 */
public abstract class ReadPairConcordanceCalculator {
	private static final Log log = Log.getInstance(ReadPairConcordanceCalculator.class);
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
					&& read1.getReferenceIndex().equals(read1.getMateReferenceIndex())
					&& read1.getReadNegativeStrandFlag() != read1.getMateNegativeStrandFlag();
		} else {
			return read1.getReadPairedFlag()
					&& !read1.getReadUnmappedFlag()
					&& !read2.getReadUnmappedFlag()
					&& read1.getReferenceIndex().equals(read2.getReferenceIndex())
					&& read1.getReadNegativeStrandFlag() != read2.getReadNegativeStrandFlag(); 
		}
	}
	public static ReadPairConcordanceCalculator create(int minFragSize, int maxFragSize, Double concordantPortion, InsertSizeDistribution insert, IdsvMetrics idsv) {
		if (maxFragSize > 0) {
			if (minFragSize > maxFragSize) {
				throw new IllegalArgumentException("minFragSize cannot be greater than maxFragSize");
			}
			return new FixedSizeReadPairConcordanceCalculator(minFragSize, maxFragSize);
		} else if (concordantPortion != null) {
			return new PercentageReadPairConcordanceCalculator(insert, concordantPortion);
		} else {
			if (idsv != null) {
				// Safety check for older versions of BWA which sets proper pair flag based only correct chromosome and orientation
				if (idsv.MAX_PROPER_PAIR_FRAGMENT_LENGTH != null &&
						idsv.MAX_PROPER_PAIR_FRAGMENT_LENGTH >= 200000 &&
						//idsv.MAX_PROPER_PAIR_FRAGMENT_LENGTH >= getMetrics().getInsertSizeMetrics().MAX_INSERT_SIZE / 2 &&
						idsv.MAX_PROPER_PAIR_FRAGMENT_LENGTH >= 10 * idsv.MAX_READ_MAPPED_LENGTH) {
					String msg = String.format("Proper pair flag indicates fragment size of %d is expected!"
									+ " This is unusually high and indicates the aligner has sets the proper pair flag based only on chromosome and orientation."
									+ " Realign with an aligner that consider fragment size when setting the proper pair flag or "
									+ " specify fixed or percentage bounds for read pair concordance.",
							idsv.MAX_PROPER_PAIR_FRAGMENT_LENGTH);
					log.error(msg);
					throw new IllegalArgumentException(msg);
				}
			}
			return new SAMFlagReadPairConcordanceCalculator(idsv);
		}
	}
}
