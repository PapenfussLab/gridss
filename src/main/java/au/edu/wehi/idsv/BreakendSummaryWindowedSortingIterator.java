package au.edu.wehi.idsv;

import au.edu.wehi.idsv.util.WindowedSortingIterator;

import java.util.Comparator;
import java.util.Iterator;

public class BreakendSummaryWindowedSortingIterator<T extends BreakendSummary> extends WindowedSortingIterator<T> {
	@SuppressWarnings("unchecked")
	public BreakendSummaryWindowedSortingIterator(final GenomicProcessingContext processContext, final int windowSize, final Iterator<T> it) {
		super(it, arg -> processContext.getLinear().getLinearCoordinate(arg.referenceIndex, arg.start), windowSize, (Comparator<T>)BreakendSummary.ByStartEnd);
	}
}
