package au.edu.wehi.idsv;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import java.nio.charset.StandardCharsets;

import org.apache.commons.lang3.StringUtils;

import au.edu.wehi.idsv.vcf.SvType;
import au.edu.wehi.idsv.vcf.VcfAttributes;
import au.edu.wehi.idsv.vcf.VcfConstants;
import au.edu.wehi.idsv.vcf.VcfSvConstants;

/**
 * Builder for generating VCF structural variation calls with appropriate attributes
 * @author Daniel Cameron
 *
 */
public class VariantContextDirectedBreakpointBuilder extends VariantContextBuilder {
	public static final String SOURCE_NAME = "idsv";
	private final ProcessingContext processContext;
	private byte[] breakendBaseQual = null;
	private void init() {
		attribute(VcfSvConstants.SV_TYPE_KEY, SvType.BND.name());
	}
	public VariantContextDirectedBreakpointBuilder(ProcessingContext processContext) {
		super();
		this.processContext = processContext;
		init();
	}
	public VariantContextDirectedBreakpointBuilder(ProcessingContext processContext, VariantContext parent) {
		super(parent);
		this.processContext = processContext;
		init();
	}
	private String getBreakendString() {
		return processContext.getVcf41Mode() ? VcfConstants.VCF41BREAKEND_REPLACEMENT : VcfConstants.VCF42BREAKEND;
	}
	/**
	 * Indicates that realignment of the breakpoint sequence failed
	 * @return builder
	 */
	public VariantContextDirectedBreakpointBuilder realignmentFailed() {
		attribute(VcfAttributes.REALIGNMENT_FAILURE.attribute(), true);
		return this;
	}
	/**
	 * Sets the variant to the given breakend.
	 * Supporting evidence is not set
	 * @param loc location of breakend
	 * @param untemplatedSequence untemplated breakend sequence 
	 * @return builder
	 */
	public VariantContextDirectedBreakpointBuilder breakend(BreakendSummary loc, String untemplatedSequence) {
		return dobreak(loc, untemplatedSequence);
	}
	/**
	 * Sets the variant to the given breakend.
	 * Supporting evidence is not set
	 * @param loc location of breakend
	 * @param untemplatedSequence untemplated breakend sequence
	 * @param  untemplatedBaseQual base qualities of untemplated breakend sequence 
	 * @return builder
	 */
	public VariantContextDirectedBreakpointBuilder breakend(BreakendSummary loc, byte[] untemplatedSequence, byte[] untemplatedBaseQual) {
		this.breakendBaseQual = untemplatedBaseQual;
		return dobreak(loc, new String(untemplatedSequence, StandardCharsets.US_ASCII));
	}
	/**
	 * Sets the variant to the given breakpoint
	 * Supporting evidence is not set
	 * @param loc location of breakpoint
	 * @param untemplatedSequence untemplated breakpoint sequence
	 * @return builder
	 */
	public VariantContextDirectedBreakpointBuilder breakpoint(BreakpointSummary loc, String untemplatedSequence) {
		return dobreak(loc, untemplatedSequence);
	}
	/**
	 * Sets the variant to the given breakend or breakpoint
	 * @param loc location of break
	 * @param untemplatedSequence untemplated break sequence
	 * @return builder
	 */
	private VariantContextDirectedBreakpointBuilder dobreak(BreakendSummary loc, String untemplatedSequence) {
		if (untemplatedSequence == null) untemplatedSequence = "";
		String chr = processContext.getDictionary().getSequence(loc.referenceIndex).getSequenceName();
		String ref, alt;
		if (processContext.getReference() == null) {
			ref = StringUtils.repeat("N", loc.end - loc.start + 1);
		} else {
			ref = new String(processContext.getReference().getSubsequenceAt(chr, loc.start, loc.start).getBases(), StandardCharsets.US_ASCII);
		}
		if (loc instanceof BreakpointSummary) {
			BreakpointSummary bp = (BreakpointSummary)loc;
			char remoteBracket = bp.direction2 == BreakendDirection.Forward ? ']' : '[';
			String target = String.format("%s:%d", processContext.getDictionary().getSequence(bp.referenceIndex2).getSequenceName(), bp.start2);
			if (loc.direction == BreakendDirection.Forward) {
				alt = String.format("%s%s%c%s%c", ref, untemplatedSequence, remoteBracket, target, remoteBracket);
			} else {
				alt = String.format("%c%s%c%s%s", remoteBracket, target, remoteBracket, untemplatedSequence, ref);
			}
		} else {
			if (loc.direction == BreakendDirection.Forward) {
				alt = ref + untemplatedSequence + getBreakendString();
			} else {
				alt = getBreakendString() + untemplatedSequence + ref;
			}
		}
		// populate
		this.loc(chr, loc.start, loc.start);
		this.alleles(ref, alt);
		this.attribute(VcfSvConstants.SV_TYPE_KEY, SvType.BND.name());
		if (loc.end != loc.start) {
			// Set confidence interval on the call if we don't have an exact breakpoint position
			this.attribute(VcfSvConstants.CONFIDENCE_INTERVAL_START_POSITION_KEY, new int[] {0, loc.end - loc.start});
		}
		if (loc instanceof BreakpointSummary) {
			BreakpointSummary bp = (BreakpointSummary)loc;
			if (bp.end2 != bp.start2) {
				this.attribute(VcfAttributes.CONFIDENCE_INTERVAL_REMOTE_BREAKEND_START_POSITION_KEY.attribute(), new int[] { 0, bp.end2 - bp.start2});
			}
		}
		return this;
	}
	/**
	 * Associates the given evidence with the variant
	 * @param e evidence to associate
	 * @return builder
	 */
	public VariantContextDirectedBreakpointBuilder evidence(EvidenceMetrics e) {
		if (e != null) {
			for (VcfAttributes a : VcfAttributes.evidenceValues()) {
				if (e.get(a) != 0) {
					// only write if we need to
					attribute(a.attribute(), e.get(a));
				} else {
					rmAttribute(a.attribute());
				}
			}
			log10PError(e.getScore() / -10);
		} else {
			log10PError(VariantContext.NO_LOG10_PERROR);
		}
		return this;
	}
	/**
	 * Associates the given assembly with the variant
	 * @param f 
	 * @return builder
	 */
	public VariantContextDirectedBreakpointBuilder assembly(
			String assemblyProgram,
			byte[] sequence,
			byte[] baseQual,
			double breakpointQuality) {
		attribute(VcfAttributes.ASSEMBLY_PROGRAM.attribute(), assemblyProgram);
		attribute(VcfAttributes.ASSEMBLY_CONSENSUS.attribute(), new String(sequence, StandardCharsets.US_ASCII));
		attribute(VcfAttributes.ASSEMBLY_QUALITY.attribute(), breakpointQuality);
		return this;
	}
	public VariantContextBuilder attribute(final VcfAttributes key, final Object value) {
		return attribute(key.attribute(), value);
	}
	@Override
	public VariantContextDirectedBreakpoint make() {
        VariantContext underlying = super.make();
        if (underlying.getEnd() < underlying.getStart() && underlying.getEnd() == -1) {
        	throw new IllegalStateException(String.format("Sanity check failure: stop not set for %s", underlying)); 
        }
		return new VariantContextDirectedBreakpoint(processContext, underlying, this.breakendBaseQual);
	}
}

