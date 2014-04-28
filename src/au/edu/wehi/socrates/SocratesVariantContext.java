package au.edu.wehi.socrates;

import net.sf.samtools.SAMSequenceDictionary;

import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;

import au.edu.wehi.socrates.vcf.VcfConstants;

/**
 * Generates variant context records from the underlying @link {@link VariantContext}
 * @author Daniel Cameron
 *
 */
public class SocratesVariantContext extends VariantContext {
	private final SAMSequenceDictionary dictionary;
	public SocratesVariantContext(SAMSequenceDictionary dictionary, VariantContext context) {
		super(context);
		this.dictionary = dictionary;
	}
	protected int getReferenceIndex(String contig) {
		return dictionary.getSequenceIndex(contig);
	}
	/**
     * @return reference index for the given sequence name, or -1 if the variant is not on a reference contig
     */
	public int getReferenceIndex() {
		return getReferenceIndex(getChr());
	}
	/**
	 * Creates a wrapper object of the appropriate type from the given {@link VariantContext} 
	 * @param context variant context
	 * @return variant context sub-type
	 */
	public static SocratesVariantContext create(SAMSequenceDictionary dictionary, VariantContext context) {
		DirectedBreakpointAssembly dba = new DirectedBreakpointAssembly(dictionary, context);
		if (dba.isValid()) return dba;
		VariantContextDirectedBreakpoint vcdp = new VariantContextDirectedBreakpoint(dictionary, context);
		if (vcdp.isValid()) return vcdp;
		// Not a variant generated or handled by us
		return new SocratesVariantContext(dictionary, context);
	}
	public boolean isValid() {
		return true;
	}
	protected Integer getOptionalIntegerAttribute(String key) {
		if (!hasAttribute(key)) return null;
		return getAttributeAsInt(key, 0);
	}
	public Integer getReferenceSpanningPairs() { return getOptionalIntegerAttribute(VcfConstants.REFERENCE_SPANNING_READ_PAIR_COUNT); }
	public Integer getReferenceReadDepth() { return getOptionalIntegerAttribute(VcfConstants.REFERENCE_READ_COUNT); }
	public Integer getOEACount() { return getOptionalIntegerAttribute(VcfConstants.UNMAPPED_MATE_READ_COUNT); }
	public Integer getSoftClipCount() { return getOptionalIntegerAttribute(VcfConstants.SOFT_CLIP_READ_COUNT); }
	public Integer getDiscordantPairCount() { return getOptionalIntegerAttribute(VcfConstants.DISCORDANT_READ_PAIR_COUNT); }
	@SuppressWarnings("unchecked")
	public <T extends SocratesVariantContext> T addEvidenceAttributes(
			T variant,
			Integer referenceSpanningPairs,
			Integer referenceReadDepth,
			Integer openEndedAnchorCount,
			Integer softClipCount,
			Integer discordantPairs) {
		VariantContextBuilder builder = new VariantContextBuilder(variant)
			.attribute(VcfConstants.REFERENCE_SPANNING_READ_PAIR_COUNT, referenceSpanningPairs)
			.attribute(VcfConstants.REFERENCE_READ_COUNT, referenceReadDepth)
			.attribute(VcfConstants.UNMAPPED_MATE_READ_COUNT, openEndedAnchorCount)
			.attribute(VcfConstants.SOFT_CLIP_READ_COUNT, softClipCount)
			.attribute(VcfConstants.DISCORDANT_READ_PAIR_COUNT, discordantPairs);
		return (T)create(variant.dictionary, builder.make());
	}
}
