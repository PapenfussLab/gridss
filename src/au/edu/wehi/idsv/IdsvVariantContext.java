package au.edu.wehi.idsv;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.google.common.base.Function;
import com.google.common.collect.ComparisonChain;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.collect.Ordering;
import com.google.common.collect.Sets;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Floats;
import com.google.common.primitives.Ints;

/**
 * Generates variant context records from the underlying @link {@link VariantContext}
 * @author Daniel Cameron
 *
 */
public class IdsvVariantContext extends VariantContext {
	protected final ProcessingContext processContext;
	protected final EvidenceSource source;
	public IdsvVariantContext(ProcessingContext processContext, EvidenceSource source, VariantContext context) {
		super(context);
		this.processContext = processContext;
		this.source = source;
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
	 * Gets the source of this evidence
	 * @return EvidenceSource if from a single source, null otherwise
	 */
	public EvidenceSource getEvidenceSource() {
		return source;
	}
	/**
	 * Creates a wrapper object of the appropriate type from the given {@link VariantContext} 
	 * @param context variant context
	 * @return variant context sub-type
	 */
	public static IdsvVariantContext create(ProcessingContext processContext, EvidenceSource source, VariantContext variant) {
		VcfBreakendSummary vbs = new VcfBreakendSummary(processContext, variant);
		if (vbs.location instanceof BreakpointSummary) return new VariantContextDirectedBreakpoint(processContext, source, variant);
		if (vbs.location instanceof BreakendSummary) return new VariantContextDirectedEvidence(processContext, source, variant);
		// Not a SV variant we're interested in
		return new IdsvVariantContext(processContext, source, variant);
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
	protected double getAttributeAsDoubleListOffset(String attrName, int offset, double defaultValue) {
		List<Double> list = getAttributeAsDoubleList(this, attrName, ImmutableList.<Double>of(), false);
		if (list == null || list.size() <= offset || list.get(offset) == null) return defaultValue;
		return list.get(offset);
	}
	protected List<Double> getAttributeAsDoubleList(String attrName) {
		List<Double> list = getAttributeAsDoubleList(this, attrName, ImmutableList.<Double>of(), true);
		return list;
	}
	protected static List<Double> getAttributeAsDoubleList(VariantContext variant, String attrName, List<Double> defaultValue, boolean throwOnParseError) {
		List<Float> floatList = getAttributeAsFloatList(variant, attrName, null, false);
		if (floatList != null) {
			List<Double> doubleList = Lists.transform(floatList, new Function<Float, Double>() {
				@Override
				public Double apply(Float arg0) {
					return (Double)(double)(float)arg0;
				}});
			return doubleList;
		}
		Object x = variant.getAttribute(attrName);
		if ( x == null || x == VCFConstants.MISSING_VALUE_v4 ) return defaultValue;
		if (x instanceof Double) {
			return Doubles.asList((Double)x);
		} else if (x instanceof double[]) {
			return Doubles.asList((double[])variant.getAttribute(attrName));
		} else if (x instanceof Iterable<?>) {
			List<Double> list = Lists.newArrayList();
			for (Object o : (Iterable<?>)x) {
				if (o instanceof Double) {
					list.add((Double)o);
				} else if (o instanceof String) {
					list.add(Double.parseDouble((String)o));
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
	private static List<Float> getAttributeAsFloatList(VariantContext variant, String attrName, List<Float> defaultValue, boolean throwOnParseError) {
		Object x = variant.getAttribute(attrName);
		if ( x == null || x == VCFConstants.MISSING_VALUE_v4 ) return defaultValue;
		if (x instanceof Float) {
			return Floats.asList((Float)x);
		} else if (x instanceof float[]) {
			return Floats.asList((float[])variant.getAttribute(attrName));
		} else if (x instanceof Iterable<?>) {
			List<Float> list = Lists.newArrayList();
			for (Object o : (Iterable<?>)x) {
				if (o instanceof Float) {
					list.add((Float)o);
				} else if (o instanceof String) {
					list.add(Float.parseFloat((String)o));
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
	protected List<String> getAttributeAsStringList(String attrName) {
		return getAttributeAsStringList(this, attrName, ImmutableList.<String>of());
	}
	protected static List<String> getAttributeAsStringList(VariantContext variant, String attrName, List<String> defaultValue) {
		Object x = variant.getAttribute(attrName);
		if ( x == null || x == VCFConstants.MISSING_VALUE_v4 ) return defaultValue;
		if (x instanceof String) {
			return ImmutableList.of((String)x);
		} else if (x instanceof String[]) {
			return Arrays.asList((String[])x);
		} else if (x instanceof Iterable<?>) {
			List<String> list = Lists.newArrayList();
			for (Object o : (Iterable<?>)x) {
				list.add(o.toString());
			}
			return list;
		} else {
			return ImmutableList.of(x.toString());
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
