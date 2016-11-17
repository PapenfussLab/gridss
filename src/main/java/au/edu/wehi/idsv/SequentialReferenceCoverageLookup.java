package au.edu.wehi.idsv;

import java.io.Closeable;
import java.util.Iterator;
import java.util.List;
import java.util.PriorityQueue;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import com.google.common.collect.PeekingIterator;
import com.google.common.collect.Queues;

import au.edu.wehi.idsv.util.SlidingWindowList;
import au.edu.wehi.idsv.visualisation.TrackedBuffer;
import gridss.analysis.IdsvMetrics;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.filter.AggregateFilter;
import htsjdk.samtools.filter.AlignedFilter;
import htsjdk.samtools.filter.DuplicateReadFilter;
import htsjdk.samtools.filter.FilteringIterator;

/**
 * Counts the number of reads and read pairs providing support for the
 * absence of a structural variation at a given position
 * 
 * @author Daniel Cameron
 *
 */
public class SequentialReferenceCoverageLookup implements Closeable, ReferenceCoverageLookup, TrackedBuffer {
	private final List<Closeable> toClose = Lists.newArrayList();
	private final PeekingIterator<SAMRecord> reads;
	private final ReadPairConcordanceCalculator pairing;
	private final PriorityQueue<Integer> currentReferenceRead = Queues.newPriorityQueue();
	private final PriorityQueue<Integer> currentStartReferencePairs = Queues.newPriorityQueue();
	private final PriorityQueue<Integer> currentEndReferencePairs = Queues.newPriorityQueue();
	/**
	 * Maximum distance from read alignment start to last concordant support position 
	 */
	private final int maxEvidenceWindow;
	private int currentReferenceIndex = -1;
	private int currentPosition;
	private int largestWindow;
	private SlidingWindowList<Integer> readCounts;
	private SlidingWindowList<Integer> pairCounts;
	/**
	 * Used to check the data is sequential
	 */
	private SAMRecord lastRead;
	/**
	 * Creates a reference lookup from the given reads
	 * @param context processing context
	 * @param windowSize window size of out-of-order querying.
	 * @param maxFragmentSize maximum fragment
	 * @param reads reads to process. <b>Must</b> be coordinate sorted
	 */
	public SequentialReferenceCoverageLookup(Iterator<SAMRecord> it, IdsvMetrics metrics, ReadPairConcordanceCalculator pairing, int windowSize) {
		this.pairing = pairing;
		if (it instanceof Closeable) toClose.add((Closeable)it);
		this.reads = Iterators.peekingIterator(new FilteringIterator(it, new AggregateFilter(ImmutableList.of(new AlignedFilter(true), new DuplicateReadFilter()))));
		this.largestWindow = windowSize;
		this.maxEvidenceWindow = Math.max(metrics.MAX_READ_LENGTH, Math.max(metrics.MAX_READ_MAPPED_LENGTH, pairing != null ? pairing.maxConcordantFragmentSize() : 0));
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
	/* (non-Javadoc)
	 * @see au.edu.wehi.idsv.ReferenceCoverageLookup#readsSupportingNoBreakendAfter(int, int)
	 */
	@Override
	public int readsSupportingNoBreakendAfter(int referenceIndex, int position) {
		ensure(referenceIndex, position);
		return getCount(readCounts, referenceIndex, position);
	}
	/* (non-Javadoc)
	 * @see au.edu.wehi.idsv.ReferenceCoverageLookup#readPairsSupportingNoBreakendAfter(int, int)
	 */
	@Override
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
		if (currentReferenceIndex == referenceIndex) {
			if (position < currentPosition - largestWindow) {
				throw new IllegalArgumentException(String.format("Unable to rewind from position %d to %d", currentPosition, position));
			}
			if (position <= currentPosition) {
				// already processed
				return;
			}
		} else {
			if (currentReferenceIndex > referenceIndex) throw new IllegalArgumentException(String.format("Unable to rewind from reference index %d to %d", currentReferenceIndex, referenceIndex));
			currentReferenceIndex = referenceIndex;
			currentPosition = 0;
			currentReferenceRead.clear();
			currentStartReferencePairs.clear();
			currentEndReferencePairs.clear();
			readCounts = new SlidingWindowList<Integer>(largestWindow);
			pairCounts = new SlidingWindowList<Integer>(largestWindow);
		}
		// skip until we're close to out window
		while (reads.hasNext() && reads.peek().getReferenceIndex() < currentReferenceIndex) {
			checkOrdered(reads.next());
		}
		while (reads.hasNext() && reads.peek().getReferenceIndex() == currentReferenceIndex && reads.peek().getAlignmentStart() < position - largestWindow - maxEvidenceWindow) {
			checkOrdered(reads.next());
		}
		// track evidence that could be in our window 
		while (reads.hasNext() && reads.peek().getReferenceIndex() == currentReferenceIndex && reads.peek().getAlignmentStart() < position - largestWindow) {
			addRead(checkOrdered(reads.next()));
		}
		// Call all position not previously called up to largestWindow bases before our target position
		for (currentPosition = Math.max(currentPosition + 1, position - largestWindow); currentPosition <= position; currentPosition++) {
			while (reads.hasNext() && reads.peek().getReferenceIndex() == currentReferenceIndex && reads.peek().getAlignmentStart() == currentPosition) {
				addRead(checkOrdered(reads.next()));
			}
			flushQueues();
			readCounts.set(currentPosition, currentReferenceRead.size());
			pairCounts.set(currentPosition, currentEndReferencePairs.size() - currentStartReferencePairs.size());
		}
		currentPosition--;
	}
	private SAMRecord checkOrdered(SAMRecord read) {
		if (lastRead != null && (read.getReferenceIndex() < lastRead.getReferenceIndex() ||
				(read.getReferenceIndex() == lastRead.getReferenceIndex() && read.getAlignmentStart() < lastRead.getAlignmentStart()))) {
			throw new IllegalStateException(String.format("Input is not sorted read %s at %s:%d before read %s at %s:%d",
					lastRead.getReadName(),
					lastRead.getReferenceName(),
					lastRead.getAlignmentStart(),
					read.getReadName(),
					read.getReferenceName(),
					read.getAlignmentStart()));
		}
		lastRead = read;
		return read;
	}
	private void addRead(SAMRecord read) {
		if (read.getReadUnmappedFlag()) return;
		// TODO: process CIGAR instead of just taking the whole alignment length as support for the reference
		currentReferenceRead.add(read.getAlignmentEnd());
		if (isLowerMappedOfNonOverlappingConcordantPair(read)) {
			currentStartReferencePairs.add(read.getAlignmentEnd());
			currentEndReferencePairs.add(read.getMateAlignmentStart());
		}
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
	private boolean isLowerMappedOfNonOverlappingConcordantPair(SAMRecord read) {
		return !read.getReadUnmappedFlag()
				&& read.getReadPairedFlag()
				&& !read.getMateUnmappedFlag()
				&& read.getAlignmentEnd() < read.getMateAlignmentStart()
				&& read.getReferenceIndex() == read.getMateReferenceIndex()
				&& (read.getAlignmentStart() < read.getMateAlignmentStart()
						|| (read.getAlignmentStart() == read.getMateAlignmentStart() && read.getFirstOfPairFlag()))
				&& pairing.isConcordant(read);
	}
	private String trackedBufferName_currentReferenceRead = "coverage.currentReferenceRead";
	private String trackedBufferName_currentStartReferencePairs = "coverage.currentStartReferencePairs";
	private String trackedBufferName_currentEndReferencePairs = "coverage.currentEndReferencePairs";
	@Override
	public void setTrackedBufferContext(String context) {
		this.trackedBufferName_currentReferenceRead = context + ".coverage.currentReferenceRead";
		this.trackedBufferName_currentStartReferencePairs = context + ".coverage.currentStartReferencePairs";
		this.trackedBufferName_currentEndReferencePairs = context + ".coverage.currentEndReferencePairs";
	}
	@Override
	public List<NamedTrackedBuffer> currentTrackedBufferSizes() {
		return ImmutableList.of(
				new NamedTrackedBuffer(trackedBufferName_currentReferenceRead, currentReferenceRead.size()),
				new NamedTrackedBuffer(trackedBufferName_currentStartReferencePairs, currentStartReferencePairs.size()),
				new NamedTrackedBuffer(trackedBufferName_currentEndReferencePairs, currentEndReferencePairs.size())
				);
	}
}
