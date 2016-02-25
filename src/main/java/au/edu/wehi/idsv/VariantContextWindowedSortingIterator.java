package au.edu.wehi.idsv;

import java.util.Comparator;
import java.util.Iterator;

import au.edu.wehi.idsv.util.WindowedSortingIterator;

import com.google.common.base.Function;

/**
 * Sorts mostly-sorted variants according to their VCF start position
 *  
 * @author Daniel Cameron
 *
 * @param <T>
 */
public class VariantContextWindowedSortingIterator<T extends IdsvVariantContext> extends WindowedSortingIterator<T> {
	@SuppressWarnings("unchecked")
	public VariantContextWindowedSortingIterator(final GenomicProcessingContext processContext, final int windowSize, final Iterator<T> it) {
		super(it, new Function<T, Long>() {
			public Long apply(T arg) {
				return processContext.getLinear().getLinearCoordinate(arg.getContig(), arg.getStart());
			}
		}, windowSize, (Comparator<T>)IdsvVariantContext.ByLocationStart);
	}
}
