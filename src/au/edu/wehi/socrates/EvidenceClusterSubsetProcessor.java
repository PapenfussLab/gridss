package au.edu.wehi.socrates;

import java.util.Iterator;

import htsjdk.samtools.util.RuntimeEOFException;
import au.edu.wehi.socrates.graph.TrapezoidGraph;
import au.edu.wehi.socrates.graph.TrapezoidGraphNode;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;
import com.google.common.collect.Ordering;
import com.google.common.collect.UnmodifiableIterator;
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
	protected boolean filterOut(BreakpointLocation loc) {
		return loc.referenceIndex != fromReferenceIndex && loc.referenceIndex != toReferenceIndex;
	}
	@Override
	protected boolean filterOut(BreakpointInterval loc) {
		return !(loc.referenceIndex == fromReferenceIndex && loc.referenceIndex2 == toReferenceIndex)
				|| (loc.referenceIndex == toReferenceIndex && loc.referenceIndex2 == fromReferenceIndex);
	}
	@Override
	public Iterator<BreakpointLocation> iterator() {
		return new BoundsAssertionIterator(super.iterator());
	}
	/**
	 * Checks that results are on the expected chromosomes
	 * @author Daniel Cameron
	 *
	 */
	private class BoundsAssertionIterator implements Iterator<BreakpointLocation> {
		private final Iterator<BreakpointLocation> it;
		public BoundsAssertionIterator(Iterator<BreakpointLocation> it) {
			this.it = it;
		}
		@Override
		public boolean hasNext() {
			return it.hasNext();
		}
		@Override
		public BreakpointLocation next() {
			BreakpointLocation loc = it.next();
			if (loc.referenceIndex != fromReferenceIndex && loc.referenceIndex != toReferenceIndex) {
				throw new RuntimeException(String.format("Sanity check failure: breakpoint %s not on referenceIndex %d or %d.", loc, fromReferenceIndex, toReferenceIndex));
			}
			if (loc instanceof BreakpointInterval) {
				BreakpointInterval interval = (BreakpointInterval)loc;
				if (interval.referenceIndex2 != fromReferenceIndex && interval.referenceIndex2 != toReferenceIndex) {
					throw new RuntimeException(String.format("Sanity check failure: breakpoint %s not on referenceIndex %d or %d.", loc, fromReferenceIndex, toReferenceIndex));
				}
			}
			return loc;
		}
		@Override
		public void remove() {
			throw new IllegalStateException();
		}
	}
}
 