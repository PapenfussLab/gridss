package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import au.edu.wehi.idsv.sam.SAMRecordUtil;

import com.google.common.collect.PeekingIterator;

/**
 * Iterates over read pairs ordered according to the start position of the locally mapped read in the pair
 * <p>Discordantly mapped read pairs will be iterated over twice, once of each read in the pair.</p>
 * @author Daniel Cameron
 *
 */
public class SequentialNonReferenceReadPairFactory extends SequentialSAMRecordFactoryBase<SAMRecord> {
	/**
	 * <p>mate iterator <b>must</b> be mate-coordinate sorted. @see SAMRecordMateCoordinateComparator<p>
	 * @param reads
	 * @param mate
	 */
	public SequentialNonReferenceReadPairFactory(
		PeekingIterator<SAMRecord> mates) {
		super(mates);
	}
	public NonReferenceReadPair createNonReferenceReadPair(SAMRecord record, SAMEvidenceSource source) {
		if (!record.getReadPairedFlag()) return null;
		SAMRecord mate = findAssociatedSAMRecord(record);
		return NonReferenceReadPair.create(record, mate, source);
	}
	@Override
	public SAMRecord findAssociatedSAMRecord(SAMRecord record) {
		if (record == null) return null;
		if (record.getReadUnmappedFlag()) return null;
		return findMatching(record.getReferenceIndex(), record.getAlignmentStart(), record.getReadName() + (record.getFirstOfPairFlag() ? SAMRecordUtil.FIRST_OF_PAIR_NAME_SUFFIX : SAMRecordUtil.SECOND_OF_PAIR_NAME_SUFFIX));
	}
	@Override
	protected String getMatchingKey(SAMRecord record) {
		return record.getReadName() + (record.getFirstOfPairFlag() ? SAMRecordUtil.SECOND_OF_PAIR_NAME_SUFFIX : SAMRecordUtil.FIRST_OF_PAIR_NAME_SUFFIX);
	}
	@Override
	protected int getReferenceIndex(SAMRecord record) {
		return record.getMateReferenceIndex();
	}
	@Override
	protected int getPosition(SAMRecord record) {
		return record.getMateAlignmentStart();
	}
}
