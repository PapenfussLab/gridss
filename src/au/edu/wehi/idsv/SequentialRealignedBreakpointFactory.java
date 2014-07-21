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
public class SequentialRealignedBreakpointFactory extends SequentialSAMRecordFactoryBase<String> {
	/**
	 * <p>read iterator <b>must</b> be coordinate sorted<p>
	 * <p>mate iterator <b>must</b> be mate-coordinate sorted. @see SAMRecordMateCoordinateComparator<p>
	 * @param reads
	 * @param mate
	 */
	public SequentialRealignedBreakpointFactory(
			PeekingIterator<SAMRecord> realigned) {
		super(realigned);
	}
	public SAMRecord findRealignedSAMRecord(DirectedEvidence source) {
		if (source == null) return null;
		return findMatching(source.getBreakendSummary().referenceIndex, source.getBreakendSummary().start, source.getEvidenceID());
	}
	@Override
	protected String getMatchingKey(SAMRecord record) {
		return BreakpointFastqEncoding.getID(record.getReadName());
	}
	@Override
	protected int getReferenceIndex(SAMRecord record) {
		return BreakpointFastqEncoding.getReferenceIndex(record.getReadName());
	}
	@Override
	protected int getPosition(SAMRecord record) {
		return BreakpointFastqEncoding.getStartPosition(record.getReadName());
	}
}
