package au.edu.wehi.idsv;

import java.util.Comparator;
import java.util.Iterator;

import au.edu.wehi.idsv.util.WindowedSortingIterator;

import com.google.common.base.Function;

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
public class DirectEvidenceWindowedSortingIterator<T extends DirectedEvidence> extends WindowedSortingIterator<T> {
	@SuppressWarnings("unchecked")
	public DirectEvidenceWindowedSortingIterator(final ProcessingContext processContext, final int windowSize, final Iterator<T> it) {
		super(it, new Function<T, Long>() {
			public Long apply(T arg) {
				return processContext.getLinear().getStartLinearCoordinate(arg.getBreakendSummary());
			}
		}, windowSize, (Comparator<T>)DirectedEvidenceOrder.ByNatural);
	}
}
