package au.edu.wehi.socrates;

import java.nio.charset.StandardCharsets;

import org.apache.commons.lang3.StringUtils;

import com.sun.javafx.collections.SetAdapterChange;

import htsjdk.samtools.SAMRecord;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import au.edu.wehi.socrates.vcf.EvidenceAttributes;
import au.edu.wehi.socrates.vcf.SvType;
import au.edu.wehi.socrates.vcf.VcfConstants;
import au.edu.wehi.socrates.vcf.VcfStructuralVariantHeaderLines;
import au.edu.wehi.socrates.vcf.VcfSvConstants;

/**
 * Builder for generating VCF structural variation calls with appropriate attributes
 * @author Daniel Cameron
 *
 */
public class VariantContextDirectedBreakpointBuilder extends VariantContextBuilder {
	public static final String SOURCE_NAME = "socrates";
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
	/**
	 * Indicates that realignment of the breakpoint sequence failed
	 * @return builder
	 */
	public VariantContextDirectedBreakpointBuilder realignmentFailed() {
		attribute(VcfConstants.REALIGNMENT_FAILURE, true);
		return this;
	}
	/**
	 * Sets the VCF attributes for a breakend at the given location
	 * @param loc
	 * @return
	 */
	private String setVcfLocationAsStartOfInterval(BreakendSummary loc) {
		this.chr(processContext.getDictionary().getSequence(loc.referenceIndex).getSequenceName());
		this.start(loc.start);
		this.stop(loc.start);
		if (loc.end != loc.start) {
			// Set confidence interval on the call if we don't have an exact breakpoint position
			this.attribute(VcfSvConstants.CONFIDENCE_INTERVAL_START_POSITION_KEY, 0);
			this.attribute(VcfSvConstants.CONFIDENCE_INTERVAL_END_POSITION_KEY, loc.end - loc.start);
		}
		String chr = processContext.getDictionary().getSequence(loc.referenceIndex).getSequenceName();
		String refBases;
		if (processContext.getReference() == null) {
			refBases = StringUtils.repeat("N", loc.end - loc.start + 1);
		} else {
			refBases = new String(processContext.getReference().getSubsequenceAt(chr, loc.start, loc.start).getBases(), StandardCharsets.US_ASCII);
		}
		return refBases; 
	}
	/**
	 * Sets the variant to the given breakend.
	 * Supporting evidence is not set
	 * @param loc location of breakend
	 * @param untemplatedSequence untemplated breakend sequence 
	 * @return builder
	 */
	public VariantContextDirectedBreakpointBuilder breakend(BreakendSummary loc, String untemplatedSequence) {
		if (loc instanceof BreakpointSummary) {
			return breakpoint((BreakpointSummary)loc, untemplatedSequence);
		}
		if (untemplatedSequence == null) untemplatedSequence = "";
		String referenceBase = setVcfLocationAsStartOfInterval(loc);
		String alt;
		if (loc.direction == BreakendDirection.Forward) {
			alt = referenceBase + untemplatedSequence + '.';
		} else {
			alt = '.' + untemplatedSequence + referenceBase;
		}
		alleles(referenceBase, alt);
		if (loc.end != loc.start) {
			attribute(VcfSvConstants.CONFIDENCE_INTERVAL_START_POSITION_KEY, new int[] { 0, loc.end - loc.start });
		}
		attribute(VcfSvConstants.SV_TYPE_KEY, SvType.BND.name());
		return this;
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
		breakend(loc, new String(untemplatedSequence, StandardCharsets.US_ASCII));
		this.breakendBaseQual = untemplatedBaseQual;
		return this;
	}
	/**
	 * Sets the variant to the given breakpoint
	 * Supporting evidence is not set
	 * @param loc location of breakpoint
	 * @param untemplatedSequence untemplated breakpoint sequence
	 * @return builder
	 */
	public VariantContextDirectedBreakpointBuilder breakpoint(BreakpointSummary loc, String untemplatedSequence) {
		if (untemplatedSequence == null) untemplatedSequence = "";
		setVcfLocationAsStartOfInterval(loc);
		String referenceBase = setVcfLocationAsStartOfInterval(loc);	
		char remoteBracket = loc.direction2 == BreakendDirection.Forward ? ']' : '[';
		String target = String.format("%s:%d", processContext.getDictionary().getSequence(loc.referenceIndex2).getSequenceName(), loc.start2);
		if (loc.direction == BreakendDirection.Forward) {
			alleles(referenceBase, String.format("%s%s%c%s%c", referenceBase, untemplatedSequence, remoteBracket, target, remoteBracket));
		} else {
			alleles(referenceBase, String.format("%c%s%c%s%s", remoteBracket, target, remoteBracket, untemplatedSequence, referenceBase));
		}
		attribute(VcfSvConstants.SV_TYPE_KEY, SvType.BND.name());
		return this;
	}
	/**
	 * Associates the given evidence with the variant
	 * @param e evidence to associate
	 * @return builder
	 */
	public VariantContextDirectedBreakpointBuilder evidence(EvidenceMetrics e) {
		if (e != null) {
			for (EvidenceAttributes a : EvidenceAttributes.values()) {
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
		attribute(VcfConstants.ASSEMBLY_PROGRAM, assemblyProgram);
		attribute(VcfConstants.ASSEMBLY_CONSENSUS, new String(sequence, StandardCharsets.US_ASCII));
		attribute(VcfConstants.ASSEMBLY_QUALITY, breakpointQuality);
		return this;
	}
	@Override
	public VariantContextDirectedBreakpoint make() {
		return new VariantContextDirectedBreakpoint(processContext, super.make(), this.breakendBaseQual);
	}
}

