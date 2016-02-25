package au.edu.wehi.idsv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import au.edu.wehi.idsv.vcf.VcfAttributes;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Ordering;
import com.google.common.primitives.Doubles;

/**
 * Generates variant context records from the underlying @link {@link VariantContext}
 * @author Daniel Cameron
 *
 */
public class IdsvVariantContext extends VariantContext {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	protected final GenomicProcessingContext processContext;
	protected final EvidenceSource source;
	public IdsvVariantContext(GenomicProcessingContext processContext, EvidenceSource source, VariantContext context) {
		super(context);
		this.processContext = processContext;
		this.source = source;
	}
	protected int getAttributeIntSum(VcfAttributes attr) { return AttributeConverter.asIntList(getAttribute(attr.attribute())).stream().mapToInt(x -> x).sum(); }
	protected int getAttributeIntOffset(VcfAttributes attr, int offset) { return AttributeConverter.asIntListOffset(getAttribute(attr.attribute()), offset, 0); }
	protected double getAttributeDoubleSum(VcfAttributes attr) { return AttributeConverter.asDoubleList(getAttribute(attr.attribute())).stream().mapToDouble(x -> x).sum(); }
	protected double getAttributeDoubleOffset(VcfAttributes attr, int offset) { return AttributeConverter.asDoubleListOffset(getAttribute(attr.attribute()), offset, 0); }
	/**
     * @return reference index for the given sequence name, or -1 if the variant is not on a reference contig
     */
	public int getReferenceIndex() {
		return getReferenceIndex(processContext, this);
	}
	/**
     * @return reference index for the given sequence name, or -1 if the variant is not on a reference contig
     */
	public static int getReferenceIndex(GenomicProcessingContext processContext, VariantContext variant) {
		return processContext.getDictionary().getSequenceIndex(variant.getContig());
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
	public static IdsvVariantContext create(GenomicProcessingContext processContext, EvidenceSource source, VariantContext variant) {
		VcfBreakendSummary vbs = new VcfBreakendSummary(processContext, variant);
		if (vbs.location instanceof BreakpointSummary) return new VariantContextDirectedBreakpoint(processContext, source, variant);
		if (vbs.location instanceof BreakendSummary) return new VariantContextDirectedEvidence(processContext, source, variant);
		// Not a SV variant we're interested in
		return new IdsvVariantContext(processContext, source, variant);
	}
	public boolean isValid() {
		return true;
	}
	public static Ordering<IdsvVariantContext> ByLocationStart = new Ordering<IdsvVariantContext>() {
		public int compare(IdsvVariantContext o1, IdsvVariantContext o2) {
			return ComparisonChain.start()
			        .compare(o1.getReferenceIndex(), o2.getReferenceIndex())
			        .compare(o1.getStart(), o2.getStart())
			        .compare(o1.getEnd(), o2.getEnd())
			        .compare(o1.getID(), o2.getID())
			        .result();
		  }
	};
	public static Ordering<IdsvVariantContext> ByLocationEnd = new Ordering<IdsvVariantContext>() {
		public int compare(IdsvVariantContext o1, IdsvVariantContext o2) {
			return ComparisonChain.start()
			        .compare(o1.getReferenceIndex(), o2.getReferenceIndex())
			        .compare(o1.getEnd(), o2.getEnd())
			        .compare(o1.getStart(), o2.getStart())
			        .compare(o1.getID(), o2.getID())
			        .result();
		  }
	};
	public static Ordering<IdsvVariantContext> ByQual = new Ordering<IdsvVariantContext>() {
		public int compare(IdsvVariantContext o1, IdsvVariantContext o2) {
			return Doubles.compare(o1.getPhredScaledQual(), o2.getPhredScaledQual());
		  }
	};
	public static Ordering<IdsvVariantContext> ByID = new Ordering<IdsvVariantContext>() {
		public int compare(IdsvVariantContext o1, IdsvVariantContext o2) {
			return ComparisonChain.start()
					.compare(o1.getID(), o2.getID())
					.result();
		  }
	};
	public static Ordering<VariantContext> VariantContextByLocationStart(final SAMSequenceDictionary dictionary) {
		return new Ordering<VariantContext>() {
			public int compare(VariantContext o1, VariantContext o2) {
				return ComparisonChain.start()
				        .compare(dictionary.getSequenceIndex(o1.getContig()), dictionary.getSequenceIndex(o2.getContig()))
				        .compare(o1.getEnd(), o2.getEnd())
				        .compare(o1.getStart(), o2.getStart())
				        .compare(o1.getID(), o2.getID())
				        .result();
			  }
		};
	}	
}
