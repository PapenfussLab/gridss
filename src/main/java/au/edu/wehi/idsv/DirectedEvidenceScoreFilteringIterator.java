package au.edu.wehi.idsv;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

import java.util.Iterator;

import com.google.common.collect.AbstractIterator;

/**
 * Only includes events higher than the given threshold score
 * @author Daniel Cameron
 *
 */
public class DirectedEvidenceScoreFilteringIterator<T extends DirectedEvidence> extends AbstractIterator<T> implements CloseableIterator<T> {
	private final Iterator<T> it;
	private final double minimumBreakendScore;
	private final double minimumBreakpointScore;
	/**
	 * Creates a new filtering iterator
	 * @param it iterator to filter
	 * @param minimumBreakendScore score at which to filter breakends
	 * @param minimumBreakpointScore score at which to filter breakpoints
	 */
	public DirectedEvidenceScoreFilteringIterator(
			Iterator<T> it,
			double minimumBreakendScore,
			double minimumBreakpointScore) {
		this.it = it;
		this.minimumBreakendScore = minimumBreakendScore;
		this.minimumBreakpointScore = minimumBreakpointScore;
	}
	@Override
	protected T computeNext() {
		while (it.hasNext()) {
			T e = it.next();
			if (e instanceof DirectedBreakpoint) {
				if (((DirectedBreakpoint) e).getBreakpointQual() >= minimumBreakpointScore) {
					return e;
				}
			} else {
				if (e.getBreakendQual() >= minimumBreakendScore) {
					return e;
				}
			}
		}
		return endOfData();
	}
	@Override
	public void close() {
		CloserUtil.close(it);
	}
}
