package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;

import com.google.common.collect.PeekingIterator;

/**
 * Finds corresponding realigned breakpoint sequence SAMRecords
 * <p>Input iterator must be ordered according to the start position of the source @See DirectedBreakpoint</p>
 * 
 * @author Daniel Cameron
 *
 */
public class SequentialRealignedBreakpointFactory extends SequentialSAMRecordFactoryBase<DirectedEvidence> {
	/**
	 * <p>read iterator <b>must</b> be coordinate sorted<p>
	 * @param reads
	 * @param mate
	 */
	public SequentialRealignedBreakpointFactory(
			PeekingIterator<SAMRecord> realigned) {
		super(realigned);
	}
	public SAMRecord findAssociatedSAMRecord(DirectedEvidence source) {
		if (source == null) return null;
		return findMatching(
				BreakpointFastqEncoding.getReferenceIndex(source),
				BreakpointFastqEncoding.getStartPosition(source),
				BreakpointFastqEncoding.getID(source));
	}
	@Override
	protected String getMatchingKey(SAMRecord record) {
		return BreakpointFastqEncoding.getEncodedID(record.getReadName());
	}
	@Override
	protected int getReferenceIndex(SAMRecord record) {
		return BreakpointFastqEncoding.getEncodedReferenceIndex(record.getReadName());
	}
	@Override
	protected int getPosition(SAMRecord record) {
		return BreakpointFastqEncoding.getEncodedStartPosition(record.getReadName());
	}
}
