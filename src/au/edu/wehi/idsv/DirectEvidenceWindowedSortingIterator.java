package au.edu.wehi.idsv;

import htsjdk.samtools.util.Log;

import java.util.Iterator;
import java.util.PriorityQueue;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

/**
 * Sorts directed evidence within from a sequence where the sequence position of
 * evidence is under a fixed distance from the sorted position
 * 
 * As SAM/BAM input is sorted by alignment start position, sorting on evidence
 * position does not require a full sort as the difference between breakend
 * start position and the alignment start position is bounded by the fragment size
 * for read pair evidence, and the read length for soft clip evidence.
 * 
 * @author Daniel Cameron
 *
 * @param <T>
 */
public class DirectEvidenceWindowedSortingIterator<T extends DirectedEvidence> extends AbstractIterator<T> {
	private static final Log log = Log.getInstance(DirectEvidenceWindowedSortingIterator.class);
	private static final int INITIAL_BUFFER_SIZE = 32;
	private final PriorityQueue<T> calls = new PriorityQueue<T>(INITIAL_BUFFER_SIZE, DirectedEvidenceOrder.ByNatural);
	private final ProcessingContext processContext;
	private final int windowSize;
	private final PeekingIterator<T> it;
	private long lastPosition = Long.MIN_VALUE;
	public DirectEvidenceWindowedSortingIterator(ProcessingContext processContext, int windowSize, Iterator<T> it) {
		this.processContext = processContext;
		this.windowSize = windowSize;
		this.it = Iterators.peekingIterator(it);
	}
	@Override
	protected T computeNext() {
		advanceUnderlying();
		if (calls.isEmpty()) return endOfData();
		T next = calls.poll();
		long nextPos = processContext.getLinear().getStartLinearCoordinate(next.getBreakendSummary());
		if (nextPos < lastPosition) {
			log.error("Sanity check failure: sorting window size too small: evidence out of order at " + next.getBreakendSummary().toString(processContext));
		}
		return next;
	}
	private void advanceUnderlying() {
		while (it.hasNext() && (calls.isEmpty() || nextRecordCouldBeAtStartOfWindow())) {
			T next = it.next();
			if (next == null) {
				throw new RuntimeException("Sanity check failure: null evidence");
			}
			calls.add(next);
		}
	}
	private boolean nextRecordCouldBeAtStartOfWindow() {
		long bufferPosition = processContext.getLinear().getStartLinearCoordinate(calls.peek().getBreakendSummary());
		long nextPosition = processContext.getLinear().getStartLinearCoordinate(it.peek().getBreakendSummary());
		return nextPosition <= bufferPosition + windowSize;
	}
}
