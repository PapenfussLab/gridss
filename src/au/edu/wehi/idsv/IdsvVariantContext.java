package au.edu.wehi.idsv;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;

import java.util.List;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.collect.Ordering;
import com.google.common.primitives.Ints;

/**
 * Generates variant context records from the underlying @link {@link VariantContext}
 * @author Daniel Cameron
 *
 */
public class IdsvVariantContext extends VariantContext {
	protected final ProcessingContext processContext;
	private final EvidenceSource source;
	public IdsvVariantContext(ProcessingContext processContext, EvidenceSource source, VariantContext context) {
		super(context);
		this.processContext = processContext;
		this.source = source;
	}
	/**
     * @return reference index for the given sequence name, or -1 if the variant is not on a reference contig
     */
	public int getReferenceIndex() {
		return processContext.getDictionary().getSequenceIndex(getChr());
	}
	/**
	 * Creates a wrapper object of the appropriate type from the given {@link VariantContext} 
	 * @param context variant context
	 * @return variant context sub-type
	 */
	public static IdsvVariantContext create(ProcessingContext processContext, EvidenceSource source, VariantContext context) {
		VariantContextDirectedBreakpoint vcdp = new VariantContextDirectedBreakpoint(processContext, source, context);
		if (vcdp.isValid()) return vcdp;
		// Not a variant generated or handled by us
		return new IdsvVariantContext(processContext, source, context);
	}
	public boolean isValid() {
		return true;
	}
	public static Ordering<? super IdsvVariantContext> ByStartStopReferenceOrder = new Ordering<IdsvVariantContext>() {
		public int compare(IdsvVariantContext o1, IdsvVariantContext o2) {
			  return ComparisonChain.start()
			        .compare(o1.getReferenceIndex(), o2.getReferenceIndex())
			        .compare(o1.getStart(), o2.getStart())
			        .compare(o1.getEnd(), o2.getEnd())
			        .result();
		  }
	};
	protected List<Integer> parseIntList(String attrName) {
		Object x = getAttribute(attrName);
		if ( x == null || x == VCFConstants.MISSING_VALUE_v4 ) return ImmutableList.of();
		if (x instanceof Integer) {
			return Ints.asList((Integer)x);
		} else if (x instanceof int[]) {
			return Ints.asList((int[])getAttribute(attrName));
		} else if (x instanceof Iterable<?>) {
			List<Integer> list = Lists.newArrayList();
			for (Object o : (Iterable<?>)x) {
				if (o instanceof Integer) {
					list.add((Integer)o);
				} else if (o instanceof String) {
					list.add(Integer.parseInt((String)o));
				} else {
					throw new IllegalStateException(String.format("Error parsing attribute %s=%s of %s", attrName, x, super.toString()));
				}
			}
			return list;
		} else {
			throw new IllegalStateException(String.format("Error parsing attribute %s=%s of %s", attrName, x, super.toString()));
		}
	}
	public EvidenceSource getEvidenceSource() {
		return source;
	}
	public static Ordering<IdsvVariantContext> ByLocation = new Ordering<IdsvVariantContext>() {
		public int compare(IdsvVariantContext o1, IdsvVariantContext o2) {
			return ComparisonChain.start()
			        .compare(o1.getReferenceIndex(), o2.getReferenceIndex())
			        .compare(o1.start, o2.start)
			        .compare(o1.stop, o2.stop)
			        .result();
		  }
	};
}
