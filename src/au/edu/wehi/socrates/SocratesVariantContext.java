package au.edu.wehi.socrates;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Ordering;

import htsjdk.variant.variantcontext.VariantContext;

/**
 * Generates variant context records from the underlying @link {@link VariantContext}
 * @author Daniel Cameron
 *
 */
public class SocratesVariantContext extends VariantContext {
	protected final ProcessingContext processContext;
	public SocratesVariantContext(ProcessingContext processContext, VariantContext context) {
		super(context);
		this.processContext = processContext;
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
	public static SocratesVariantContext create(ProcessingContext processContext, VariantContext context) {
		VariantContextDirectedBreakpoint vcdp = new VariantContextDirectedBreakpoint(processContext, context);
		if (vcdp.isValid()) return vcdp;
		// Not a variant generated or handled by us
		return new SocratesVariantContext(processContext, context);
	}
	public boolean isValid() {
		return true;
	}
	public static Ordering<? super SocratesVariantContext> ByStartStopReferenceOrder = new Ordering<SocratesVariantContext>() {
		public int compare(SocratesVariantContext o1, SocratesVariantContext o2) {
			  return ComparisonChain.start()
			        .compare(o1.getReferenceIndex(), o2.getReferenceIndex())
			        .compare(o1.getStart(), o2.getStart())
			        .compare(o1.getEnd(), o2.getEnd())
			        .result();
		  }
	};
}
