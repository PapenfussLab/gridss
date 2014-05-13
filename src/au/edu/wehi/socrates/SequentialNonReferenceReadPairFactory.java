package au.edu.wehi.socrates;

import htsjdk.samtools.SAMRecord;

import com.google.common.collect.PeekingIterator;

/**
 * Iterates over read pairs ordered according to the start position of the locally mapped read in the pair
 * <p>Discordantly mapped read pairs will be iterated over twice, once of each read in the pair.</p>
 * @author Daniel Cameron
 *
 */
public class SequentialNonReferenceReadPairFactory extends SequentialSAMRecordFactoryBase<String> {
	/**
	 * <p>mate iterator <b>must</b> be mate-coordinate sorted. @see SAMRecordMateCoordinateComparator<p>
	 * @param reads
	 * @param mate
	 */
	public SequentialNonReferenceReadPairFactory(
		PeekingIterator<SAMRecord> mates) {
		super(mates);
	}
	public NonReferenceReadPair createNonReferenceReadPair(SAMRecord record, int maxFragmentSize) {
		if (!record.getReadPairedFlag()) return null;
		if (record.getProperPairFlag()) return null;
		SAMRecord mate = findMate(record);
		if (record != null && mate != null) {
			return new NonReferenceReadPair(record, mate, maxFragmentSize); 
		}
		return null;
	}
	private SAMRecord findMate(SAMRecord record) {
		if (record == null) return null;
		if (record.getReadUnmappedFlag()) return null;
		return findMatching(record.getReferenceIndex(), record.getAlignmentStart(), record.getReadName() + (record.getFirstOfPairFlag() ? "/1" : "/2"));
	}
	@Override
	protected String getMatchingKey(SAMRecord record) {
		return record.getReadName() + (record.getFirstOfPairFlag() ? "/2" : "/1");
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
