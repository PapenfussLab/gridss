package au.edu.wehi.idsv.vcf;

import java.util.Iterator;

import com.google.common.collect.AbstractIterator;

import htsjdk.variant.variantcontext.VariantContext;

/**
 * Removes all variants failing any filter from the sequence
 * @author Daniel Cameron
 *
 * @param <T>
 */
public class PassFiltersIterator<T extends VariantContext> extends AbstractIterator<T> {
	private Iterator<T> it;
	public PassFiltersIterator(Iterator<T> it) {
		this.it = it;
	}
	@Override
	protected T computeNext() {
		while (it.hasNext()) {
			T v = it.next();
			if (!v.isFiltered()) return v;
		}
		return endOfData();
	}

}
