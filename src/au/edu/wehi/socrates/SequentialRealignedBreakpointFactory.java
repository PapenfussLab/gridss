package au.edu.wehi.socrates;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;

import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

/**
 * Finds corresponding realigned breakpoint sequence SAMRecords
 * <p>Input iterator must be ordered according to the start position of the source @See DirectedBreakpoint</p>
 * 
 * @author Daniel Cameron
 *
 */
public class SequentialRealignedBreakpointFactory {
	private final SAMSequenceDictionary dictionary;
	private final PeekingIterator<SAMRecord> realigned;
	private int currentReferenceIndex = -1; 
	private long currentPosition = -1;
	private Map<String, SAMRecord> currentReads = new HashMap<String, SAMRecord>();
	/**
	 * <p>read iterator <b>must</b> be coordinate sorted<p>
	 * <p>mate iterator <b>must</b> be mate-coordinate sorted. @see SAMRecordMateCoordinateComparator<p>
	 * @param reads
	 * @param mate
	 */
	public SequentialRealignedBreakpointFactory(
		SAMSequenceDictionary dictionary,
		Iterator<SAMRecord> realigned) {
		this.dictionary = dictionary;
		this.realigned = Iterators.peekingIterator(realigned);
	}
	public SAMRecord findRealignedSAMRecord(DirectedBreakpoint source) {
		if (source == null) return null;
		return findRealignedSAMRecord(source.getChr(), source.getStart(), source.getBreakpointID());
	}
	private SAMRecord findRealignedSAMRecord(String chr, int position, String id) {
		load(dictionary.getSequenceIndex(chr), position);
		if (currentReads.containsKey(id)) {
			return currentReads.get(id);
		}
		return null;
	}
	private void load(int referenceIndex, int alignmentStart) {
		if (referenceIndex == currentReferenceIndex && currentPosition > alignmentStart) {
			throw new RuntimeException(String.format("Sanity check failure: cannot rewind SequentialRealignedBreakpointFactory source iterator from %d:%d to %d:%d", currentReferenceIndex, currentPosition, referenceIndex, alignmentStart));
		}
		if (referenceIndex != currentReferenceIndex || currentPosition != alignmentStart) {
			currentReads.clear();
			currentReferenceIndex = referenceIndex;
			currentPosition = alignmentStart;
			while (realigned.hasNext()) {
				SAMRecord r = realigned.peek();
				String sourceID = BreakpointFastqIDEncoding.getID(r.getReadName());
				int sourceReferenceIndex = dictionary.getSequenceIndex(BreakpointFastqIDEncoding.getChr(r.getReadName()));
				long sourcePosition = BreakpointFastqIDEncoding.getPosition(r.getReadName());
				if (sourceReferenceIndex < currentReferenceIndex) {
					// skip reads that are before our position
					realigned.next();
					continue;
				} else if (sourceReferenceIndex == currentReferenceIndex) {
					if (sourcePosition < currentPosition) {
						// skip reads that are before our position
						realigned.next();
						continue;
					} else if (sourcePosition == currentPosition) {
						// add reads at current position
						currentReads.put(sourceID, realigned.next());
						continue;
					} else {
						// stop when we hit a read after our position
						break;
					}
				} else {
					// stop when we hit a read after our position
					break;
				}
			}
		}
	}
}
