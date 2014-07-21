package au.edu.wehi.idsv;

import java.util.Iterator;
/**
 * Calls breakpoints from the given evidence
 * 
 * @author Daniel Cameron
 */
public class EvidenceClusterSubsetProcessor extends EvidenceClusterProcessor {
	private int fromReferenceIndex;
	private int toReferenceIndex;
	public EvidenceClusterSubsetProcessor(
			ProcessingContext context,
			int fromReferenceIndex,
			int toReferenceIndex) {
		super(context);
		this.fromReferenceIndex = Math.min(fromReferenceIndex, toReferenceIndex);
		this.toReferenceIndex = Math.max(fromReferenceIndex, toReferenceIndex);
	}
	@Override
	protected boolean filterOut(BreakendSummary loc) {
		return loc.referenceIndex != fromReferenceIndex && loc.referenceIndex != toReferenceIndex;
	}
	@Override
	protected boolean filterOut(BreakpointSummary loc) {
		return !(loc.referenceIndex == fromReferenceIndex && loc.referenceIndex2 == toReferenceIndex)
				|| (loc.referenceIndex == toReferenceIndex && loc.referenceIndex2 == fromReferenceIndex);
	}
	@Override
	public Iterator<VariantContextDirectedEvidence> iterator() {
		return new BoundsAssertionIterator<VariantContextDirectedEvidence>(super.iterator());
	}
	/**
	 * Checks that results are on the expected chromosomes
	 * @author Daniel Cameron
	 *
	 */
	private class BoundsAssertionIterator<T extends VariantContextDirectedEvidence> implements Iterator<T> {
		private final Iterator<T> it;
		public BoundsAssertionIterator(Iterator<T> it) {
			this.it = it;
		}
		@Override
		public boolean hasNext() {
			return it.hasNext();
		}
		@Override
		public T next() {
			T value = it.next(); 
			BreakendSummary loc = value.getBreakendSummary();
			if (loc.referenceIndex != fromReferenceIndex && loc.referenceIndex != toReferenceIndex) {
				throw new RuntimeException(String.format("Sanity check failure: breakpoint %s not on referenceIndex %d or %d.", loc, fromReferenceIndex, toReferenceIndex));
			}
			if (loc instanceof BreakpointSummary) {
				BreakpointSummary interval = (BreakpointSummary)loc;
				if (interval.referenceIndex2 != fromReferenceIndex && interval.referenceIndex2 != toReferenceIndex) {
					throw new RuntimeException(String.format("Sanity check failure: breakpoint %s not on referenceIndex %d or %d.", loc, fromReferenceIndex, toReferenceIndex));
				}
			}
			return value;
		}
		@Override
		public void remove() {
			throw new IllegalStateException();
		}
	}
}
 