package au.edu.wehi.socrates;

import java.util.ArrayDeque;
import java.util.PriorityQueue;
import java.util.Queue;

import au.edu.wehi.socrates.util.NonReferenceReadPairLocalComparator;
import au.edu.wehi.socrates.util.SAMRecordEndCoordinateComparator;
import au.edu.wehi.socrates.util.SAMRecordSummary;

import com.google.common.collect.PeekingIterator;
import com.sun.jmx.remote.internal.ArrayQueue;

import net.sf.picard.util.Log;
import net.sf.samtools.SAMRecord;

/**
 * Iterator-like object that sequentially traverses a BAM file
 * providing structural variation evidence in sequential order
 * @author Daniel Cameron
 *
 */
public class SequentialSAMRecordBuffer {
	private static final int INITIAL_BUFFER_SIZE = 1024;
	private final PeekingIterator<SAMRecord> reads;
	private final SequentialNonReferenceReadPairFactory mateFactory;
	private final int referenceIndex;
	private final Queue<SAMRecord> startSC = new ArrayDeque<SAMRecord>(INITIAL_BUFFER_SIZE);
	private final Queue<NonReferenceReadPair> startRP = new ArrayDeque<NonReferenceReadPair>(INITIAL_BUFFER_SIZE);
	private final PriorityQueue<SAMRecord> endSC = new PriorityQueue<SAMRecord>(INITIAL_BUFFER_SIZE, new SAMRecordEndCoordinateComparator());
	private final PriorityQueue<NonReferenceReadPair> endRP = new PriorityQueue<NonReferenceReadPair>(INITIAL_BUFFER_SIZE, new NonReferenceReadPairLocalComparator(new SAMRecordEndCoordinateComparator()));
	private final Log log = Log.getInstance(SequentialSAMRecordBuffer.class);
	public SequentialSAMRecordBuffer(int referenceIndex, PeekingIterator<SAMRecord> reads, SequentialNonReferenceReadPairFactory mateFactory) {
		this.referenceIndex = referenceIndex;
		this.reads = reads;
		this.mateFactory = mateFactory;
		advanceToReferenceIndex();
	}
	private <T> T getNext(Queue<T> queue, int maxPosition) {
		if (queue.isEmpty()) processTo(maxPosition);
		if (!queue.isEmpty()) return queue.poll();
		return null;
	}
	/**
	 * Returns the next read soft clipped at alignment end,
	 * with read ending at or before position 
	 * @param position max ending position to return soft clipped read
	 * @return @See SAMRecord read, null if no matching read exists
	 */
	public SAMRecord getNextEndingSoftClippedRead(int position) {
		return getNext(endSC, position);
	}
	/**
	 * Returns the next read pair ending at or before the given end position
	 * whose mate was expected to be concordantly mapped after the read but
	 * was not.
	 * @param position max ending position to return read pair
	 * @return @See NonReferenceReadPair read pair, null if no matching read exists
	 */
	public NonReferenceReadPair getNextNonReferenceForwardReadPair(int position) {
		return getNext(endRP, position);
	}
	/**
	 * Returns the next read soft clipped at alignment start at or before position 
	 * @param position max position to return soft clipped read
	 * @return @See SAMRecord read, null if no matching read exists
	 */
	public SAMRecord getNextStartingSoftClippedRead(int position) {
		return getNext(startSC, position);
	}
	/**
	 * Returns the next read pair starting at or before the given position
	 * whose mate was expected to be concordantly mapped before the read but
	 * was not.
	 * @param position max position to return read pair
	 * @return @See NonReferenceReadPair read pair, null if no matching read exists
	 */
	public NonReferenceReadPair getNextNonReferenceBackwardReadPair(int position) {
		return getNext(startRP, position);
	}
	public boolean hasReads() {
		return !startSC.isEmpty() ||
				!startRP.isEmpty() ||
				!endSC.isEmpty() ||
				!endRP.isEmpty() ||
				(reads.hasNext() && reads.peek().getReferenceIndex() == referenceIndex);
	}
	private void advanceToReferenceIndex() {
		// Advance iterator to our reference index
		if (reads.hasNext() && reads.peek().getReferenceIndex() < referenceIndex) {
			int i = 0;
			while (reads.hasNext() && reads.peek().getReferenceIndex() < referenceIndex) {
				reads.next();
				i++;
			}
			log.warn(String.format("Skipped %d records when processing referenceIndex %d", i, referenceIndex));
		}
	}
	/**
	 * Process all reads up to the given position.
	 * Since the input is sorted by starting position, we'll process reads starting
	 * at that position, which also guarantees we'll have processed all reads ending
	 * at the given position (since they have to start at or before the position)
	 * @param position position to process reads until
	 */
	private void processTo(int position) {
		// TODO: what to do about unaligned reads?
		while (reads.hasNext() && reads.peek().getReferenceIndex() == referenceIndex && reads.peek().getAlignmentStart() <= position) {
			SAMRecord next = reads.next();
			processRead(next);
		}
	}
	private void processRead(SAMRecord record) {
		if (SAMRecordSummary.getStartSoftClipLength(record) > 0) {
			startSC.add(record);
		}
		if (SAMRecordSummary.getEndSoftClipLength(record) > 0) {
			endSC.add(record);
		}
		if (SAMRecordSummary.isPartOfNonReferenceReadPair(record)) {
			NonReferenceReadPair pair = mateFactory.createNonReferenceReadPair(record);
			if (pair.getBreakpointDirection() == BreakpointDirection.Forward) { 
				endRP.add(pair);
			} else {
				startRP.add(pair);
			}
		}
	}
}
