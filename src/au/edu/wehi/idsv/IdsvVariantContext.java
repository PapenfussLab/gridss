package au.edu.wehi.idsv;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.collect.Ordering;
import com.google.common.collect.Sets;
import com.google.common.primitives.Ints;

/**
 * Generates variant context records from the underlying @link {@link VariantContext}
 * @author Daniel Cameron
 *
 */
public class IdsvVariantContext extends VariantContext {
	protected final ProcessingContext processContext;
	protected Set<EvidenceSource> sourceSet = Sets.newHashSet();
	public IdsvVariantContext(ProcessingContext processContext, Set<EvidenceSource> sourceSet, VariantContext context) {
		super(context);
		this.processContext = processContext;
		this.sourceSet = sourceSet;
	}
	/**
     * @return reference index for the given sequence name, or -1 if the variant is not on a reference contig
     */
	public int getReferenceIndex() {
		return getReferenceIndex(processContext, this);
	}
	/**
     * @return reference index for the given sequence name, or -1 if the variant is not on a reference contig
     */
	public static int getReferenceIndex(ProcessingContext processContext, VariantContext variant) {
		return processContext.getDictionary().getSequenceIndex(variant.getChr());
	}
	/**
	 * Creates a wrapper object of the appropriate type from the given {@link VariantContext} 
	 * @param context variant context
	 * @return variant context sub-type
	 */
	public static IdsvVariantContext create(ProcessingContext processContext, EvidenceSource source, VariantContext variant) {
		VcfBreakendSummary vbs = new VcfBreakendSummary(processContext, variant);
		HashSet<EvidenceSource> sourceSet = Sets.newHashSet();
		sourceSet.add(source);
		if (vbs.location instanceof DirectedBreakpoint) return new VariantContextDirectedBreakpoint(processContext, sourceSet, variant, null);
		if (vbs.location instanceof DirectedEvidence) return new VariantContextDirectedEvidence(processContext, sourceSet, variant);
		// Not a SV variant we're interested in
		return new IdsvVariantContext(processContext, sourceSet, variant);
	}
	public boolean isValid() {
		return true;
	}
	protected int getAttributeAsIntListOffset(String attrName, int offset, int defaultValue) {
		List<Integer> list = getAttributeAsIntList(this, attrName, ImmutableList.<Integer>of(), false);
		if (list == null || list.size() <= offset || list.get(offset) == null) return defaultValue;
		return list.get(offset);
	}
	protected List<Integer> getAttributeAsIntList(String attrName) {
		List<Integer> list = getAttributeAsIntList(this, attrName, ImmutableList.<Integer>of(), true);
		return list;
	}
	protected static List<Integer> getAttributeAsIntList(VariantContext variant, String attrName, List<Integer> defaultValue, boolean throwOnParseError) {
		Object x = variant.getAttribute(attrName);
		if ( x == null || x == VCFConstants.MISSING_VALUE_v4 ) return defaultValue;
		if (x instanceof Integer) {
			return Ints.asList((Integer)x);
		} else if (x instanceof int[]) {
			return Ints.asList((int[])variant.getAttribute(attrName));
		} else if (x instanceof Iterable<?>) {
			List<Integer> list = Lists.newArrayList();
			for (Object o : (Iterable<?>)x) {
				if (o instanceof Integer) {
					list.add((Integer)o);
				} else if (o instanceof String) {
					list.add(Integer.parseInt((String)o));
				} else {
					if (throwOnParseError) throw new IllegalStateException(String.format("Error parsing attribute %s=%s of %s", attrName, x, variant.toString()));
					return defaultValue;
				}
			}
			return list;
		} else {
			if (throwOnParseError) throw new IllegalStateException(String.format("Error parsing attribute %s=%s of %s", attrName, x, variant.toString()));
			return defaultValue;
		}
	}
	public static Ordering<IdsvVariantContext> ByLocationStart = new Ordering<IdsvVariantContext>() {
		public int compare(IdsvVariantContext o1, IdsvVariantContext o2) {
			return ComparisonChain.start()
			        .compare(o1.getReferenceIndex(), o2.getReferenceIndex())
			        .compare(o1.start, o2.start)
			        .compare(o1.stop, o2.stop)
			        .result();
		  }
	};
	public static Ordering<IdsvVariantContext> ByLocationEnd = new Ordering<IdsvVariantContext>() {
		public int compare(IdsvVariantContext o1, IdsvVariantContext o2) {
			return ComparisonChain.start()
			        .compare(o1.getReferenceIndex(), o2.getReferenceIndex())
			        .compare(o1.stop, o2.stop)
			        .compare(o1.start, o2.start)
			        .result();
		  }
	};
}
