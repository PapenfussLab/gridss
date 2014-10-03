package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import au.edu.wehi.idsv.sam.SamTags;

import com.google.common.collect.PeekingIterator;

/**
 * Finds the corresponding soft clip SAMRecord for a realigned breakpoint
 * <p>Input iterator must be ordered according to the start position of the corresponding realigned breakpoint
 * 
 * @author Daniel Cameron
 *
 */
public class SequentialSoftClipRealignedRemoteBreakpointFactory extends SequentialSAMRecordFactoryBase<SAMRecord> {
	private final BreakendDirection direction;
	/**
	 * <p>read iterator <b>must</b> be coordinate sorted<p>
	 * @param reads
	 * @param mate
	 */
	public SequentialSoftClipRealignedRemoteBreakpointFactory(
			PeekingIterator<SAMRecord> softclipsSortedByRealignmentLocation,
			BreakendDirection direction) {
		super(softclipsSortedByRealignmentLocation);
		this.direction = direction;
	}
	public SAMRecord findAssociatedSAMRecord(SAMRecord realignedClippedBases) {
		if (realignedClippedBases == null) return null;
		return findMatching(
				realignedClippedBases.getReferenceIndex(),
				realignedClippedBases.getAlignmentStart(),
				BreakpointFastqEncoding.getEncodedID(realignedClippedBases.getReadName()));
	}
	@Override
	protected String getMatchingKey(SAMRecord record) {
		return SoftClipEvidence.getEvidenceID(direction, record);
	}
	@Override
	protected int getReferenceIndex(SAMRecord record) {
		return record.getIntegerAttribute(SamTags.REALIGNMENT_REFERENCE_INDEX);
	}
	@Override
	protected int getPosition(SAMRecord record) {
		return record.getIntegerAttribute(SamTags.REALIGNMENT_POSITION);
	}
}
