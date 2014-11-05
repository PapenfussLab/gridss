package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.filter.AggregateFilter;
import htsjdk.samtools.filter.AlignedFilter;
import htsjdk.samtools.filter.DuplicateReadFilter;
import htsjdk.samtools.filter.FilteringIterator;

import java.io.Closeable;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.PriorityQueue;

import au.edu.wehi.idsv.util.SlidingWindowList;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;
import com.google.common.collect.Queues;

/**
 * Counts the number of reads and read pairs providing support for the
 * absence of a structural variation at a given position
 * 
 * @author Daniel Cameron
 *
 */
public class SequentialReferenceCoverageLookup implements Closeable {
	private final List<Closeable> toClose = new ArrayList<>();
	private final PeekingIterator<SAMRecord> reads;
	private final PriorityQueue<Integer> currentReferenceRead = Queues.newPriorityQueue();
	private final PriorityQueue<Integer> currentStartReferencePairs = Queues.newPriorityQueue();;
	private final PriorityQueue<Integer> currentEndReferencePairs = Queues.newPriorityQueue();;
	private int currentReferenceIndex;
	private int currentPosition;
	private int largestWindow;
	private SlidingWindowList<Integer> readCounts;
	private SlidingWindowList<Integer> pairCounts;
	/**
	 * Creates a reference lookup from the given reads
	 * @param context processing context
	 * @param windowSize window size of out-of-order querying.
	 * @param reads reads to process. <b>Must</b> be coordinate sorted
	 */
	public SequentialReferenceCoverageLookup(Iterator<SAMRecord> it, int windowSize) {
		if (it instanceof Closeable) toClose.add((Closeable)it);
		this.reads = Iterators.peekingIterator(new FilteringIterator(it, new AggregateFilter(ImmutableList.of(new AlignedFilter(true), new DuplicateReadFilter()))));
		this.largestWindow = windowSize;
		jumpTo(0);
	}
	public void close() {
		for (Closeable c : toClose) {
			try {
				c.close();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		toClose.clear();
	}
	private int getCount(SlidingWindowList<Integer> counts, int referenceIndex, int position) {
		if (counts.size() <= position) return 0;
		// 10 10 0 good
		// 2 1 1 good
		// 0 1 1 bad
		if (position < counts.size() - counts.getWindowSize()) throw new IllegalArgumentException(String.format("position %d outside of window of size %d ending at position %d", position, counts.getWindowSize(), counts.size()));
		Integer count = counts.get(position);
		if (count == null) return 0;
		return (int)count;
	}
	/**
	 * Number of reference reads providing evidence against a breakend immediately after the given base
	 * @param referenceIndex contig
	 * @param position position immediate before putative breakend
	 * @return number of reads spanning the putative breakend
	 */
	public int readsSupportingNoBreakendAfter(int referenceIndex, int position) {
		ensure(referenceIndex, position);
		return getCount(readCounts, referenceIndex, position);
	}
	/**
	 * Number of read pairs providing evidence against a breakend immediately after the given base
	 * @param referenceIndex contig
	 * @param position position immediate before putative breakend
	 * @return number of read pairs spanning the putative breakend
	 */
	public int readPairsSupportingNoBreakendAfter(int referenceIndex, int position) {
		ensure(referenceIndex, position);
		return getCount(pairCounts, referenceIndex, position);
	}
	/**
	 * Ensures the given position has been processed
	 * @param referenceIndex
	 * @param position
	 */
	private void ensure(int referenceIndex, int position) {
		if (currentReferenceIndex > referenceIndex) throw new IllegalArgumentException(String.format("Unable to rewind from reference index %d to %d", currentReferenceIndex, referenceIndex));
		if (referenceIndex > currentReferenceIndex) jumpTo(referenceIndex);
		while (currentPosition < position && !doneWithContig()) { // continue while we haven't reached our target position
			// Advance and call the next position
			currentPosition++;
			flushQueues();
			processCurrentReads();
			readCounts.set(currentPosition, currentReferenceRead.size());
			pairCounts.set(currentPosition, currentEndReferencePairs.size() - currentStartReferencePairs.size());
			if (queuesAreEmpty()) {
				// short-cut advance to next read as all intervening values are going to be zero
				currentPosition = nextReadPosition() - 1;
			}
		}
	}
	private boolean doneWithContig() {
		return !haveUnprocessReads() && queuesAreEmpty();
	}
	private int nextReadPosition() {
		if (!haveUnprocessReads()) return Integer.MAX_VALUE;
		return reads.peek().getAlignmentStart();
	}
	private boolean queuesAreEmpty() {
		return currentReferenceRead.isEmpty() && currentStartReferencePairs.isEmpty() && currentEndReferencePairs.isEmpty();
	}
	private boolean haveUnprocessReads() {
		return reads.hasNext() && reads.peek().getReferenceIndex() == currentReferenceIndex;
	}
	/**
	 * Flushes the queues of all reads that no longer support the reference
	 * at the given current position
	 */
	private void flushQueues() {
		while (!currentReferenceRead.isEmpty() && currentReferenceRead.peek() <= currentPosition) currentReferenceRead.poll();
		while (!currentStartReferencePairs.isEmpty() && currentStartReferencePairs.peek() <= currentPosition) currentStartReferencePairs.poll();
		while (!currentEndReferencePairs.isEmpty() && currentEndReferencePairs.peek() <= currentPosition) currentEndReferencePairs.poll();
	}
	/**
	 * Processes all reads starting at the current position
	 */
	private void processCurrentReads() {
		while (reads.hasNext() && reads.peek().getReferenceIndex() == currentReferenceIndex && reads.peek().getAlignmentStart() == currentPosition) {
			SAMRecord read = reads.next();
			// TODO: process CIGAR instead of just taking the whole alignment length as support for the reference
			currentReferenceRead.add(read.getAlignmentEnd());
			if (read.getReadPairedFlag()
					&& !read.getMateUnmappedFlag()
					&& read.getProperPairFlag()
					&& read.getMateReferenceIndex() == currentReferenceIndex
					// only get first of the pair
					&& (read.getAlignmentStart() < read.getMateAlignmentStart()
						|| (read.getAlignmentStart() == read.getMateAlignmentStart() && read.getFirstOfPairFlag()))) {
				// we're a proper pair
				currentStartReferencePairs.add(read.getAlignmentEnd());
				currentEndReferencePairs.add(read.getMateAlignmentStart());
			}
		}
	}
	/**
	 * Stops processing the current reference index and advances to the given index
	 * @param referenceIndex
	 */
	private void jumpTo(int referenceIndex) {
		currentReferenceIndex = referenceIndex;
		currentPosition = 0;
		currentReferenceRead.clear();
		currentStartReferencePairs.clear();
		currentEndReferencePairs.clear();
		readCounts = new SlidingWindowList<Integer>(largestWindow);
		pairCounts = new SlidingWindowList<Integer>(largestWindow);
		while (reads.hasNext() && reads.peek().getReferenceIndex() < currentReferenceIndex) reads.next();
	}
}
