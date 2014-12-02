package au.edu.wehi.idsv;

import java.util.Iterator;

/**
 * Checks that results are on the expected chromosomes
 * @author Daniel Cameron
 *
 */
class DirectedEvidenceBoundsAssertionIterator<T extends DirectedEvidence> implements Iterator<T> {
	private final Iterator<T> it;
	private final int referenceIndex;
	private final int referenceIndex2;
	public DirectedEvidenceBoundsAssertionIterator(Iterator<T> it, int referenceIndex, int referenceIndex2) {
		this.it = it;
		this.referenceIndex = referenceIndex;
		this.referenceIndex2 = referenceIndex2;
	}
	@Override
	public boolean hasNext() {
		return it.hasNext();
	}
	@Override
	public T next() {
		T value = it.next(); 
		BreakendSummary loc = value.getBreakendSummary();
		if (loc.referenceIndex != referenceIndex && loc.referenceIndex != referenceIndex2) {
			throw new RuntimeException(String.format("Sanity check failure: breakpoint %s not on referenceIndex %d or %d.", loc, referenceIndex, referenceIndex2));
		}
		if (loc instanceof BreakpointSummary) {
			BreakpointSummary interval = (BreakpointSummary)loc;
			if (interval.referenceIndex2 != referenceIndex && interval.referenceIndex2 != referenceIndex2) {
				throw new RuntimeException(String.format("Sanity check failure: breakpoint %s not on referenceIndex %d or %d.", loc, referenceIndex, referenceIndex2));
			}
		}
		return value;
	}
	@Override
	public void remove() {
		throw new IllegalStateException();
	}
}
