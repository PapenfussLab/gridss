package au.edu.wehi.socrates;

import java.nio.charset.StandardCharsets;

import org.apache.commons.lang3.StringUtils;

import htsjdk.samtools.SAMRecord;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

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
	 * Updates the given breakpoint with the given realignment of the breakpoint sequence
	 * @param evidence
	 * @param realignment
	 * @return
	 */
	public VariantContextDirectedBreakpointBuilder realigned(BreakendSummary localAnchor, SAMRecord realignment) {		
		if (realignment.getReadUnmappedFlag()) {
			attribute(VcfConstants.REALIGNMENT_FAILURE, true);
		} else {
			BreakendSummary remoteAnchor;
			int untemplatedSequenceLength;
			if ((localAnchor.direction == BreakendDirection.Forward && realignment.getReadNegativeStrandFlag()) ||
					(localAnchor.direction == BreakendDirection.Backward && !realignment.getReadNegativeStrandFlag())) {
				// negative strand template match means we flip the expected direction
				// since breakend sequence is on the +ve strand,
				// realigned forward breakends on +ve stand would indicate backward breakend
				// realigned backward breakends on +ve stand would indicate forward breakend
				remoteAnchor = new BreakendSummary(realignment.getReferenceIndex(), BreakendDirection.Forward,
						realignment.getAlignmentEnd(), realignment.getAlignmentEnd(), null);
				untemplatedSequenceLength = SAMRecordUtil.getEndSoftClipLength(realignment);
			} else {
				// ACGT. => CGT breakpoint, -> ACGT.---.CGT
				remoteAnchor = new BreakendSummary(realignment.getReferenceIndex(), BreakendDirection.Backward,
						realignment.getAlignmentStart(), realignment.getAlignmentStart(), null);
				untemplatedSequenceLength = SAMRecordUtil.getStartSoftClipLength(realignment);
			}			
			String untemplatedSequence = new String(evidence.getBreakpointSequence(), StandardCharsets.US_ASCII);
			if (localAnchor.direction == BreakendDirection.Forward) {
				untemplatedSequence = untemplatedSequence.substring(0, untemplatedSequenceLength);
			} else {
				untemplatedSequence = untemplatedSequence.substring(untemplatedSequence.length() - untemplatedSequenceLength);
			}
			location(new BreakpointSummary(localAnchor, remoteAnchor, localAnchor.qual), untemplatedSequence);
		}
		return this;
	}
	public VariantContextDirectedBreakpointBuilder breakpoint(DirectedEvidence evidence) {
		if (evidence instanceof DirectedBreakpoint) {
			DirectedBreakpoint bp = (DirectedBreakpoint)evidence;
			byte[] seq = bp.getBreakpointSequence();
			byte[] qual = bp.getBreakpointQuality();
			location(evidence.getBreakendSummary(), new String(seq, StandardCharsets.US_ASCII));
		} else {
			location(evidence.getBreakendSummary(), "");
		}
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
			// TODO: correct setting of CI_POS for all directions
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
	public VariantContextDirectedBreakpointBuilder location(BreakendSummary loc, String untemplatedSequence) {
		if (loc instanceof BreakpointSummary) {
			return location((BreakpointSummary)loc, "");
		}
		String referenceBase = setVcfLocationAsStartOfInterval(loc);
		String alt;
		if (loc.direction == BreakendDirection.Forward) {
			alt = referenceBase + untemplatedSequence + '.';
		} else {
			alt = '.' + untemplatedSequence + referenceBase;
		}
		alleles(referenceBase, alt);
		log10PError(loc.qual / -10);
		if (loc.end != loc.start) {
			attribute(VcfSvConstants.CONFIDENCE_INTERVAL_START_POSITION_KEY, new int[] { 0, loc.end - loc.start });
		}
		return this;
	}
	public VariantContextDirectedBreakpointBuilder location(BreakpointSummary loc, String untemplatedSequence) {
		setVcfLocationAsStartOfInterval(loc);
		String referenceBase = setVcfLocationAsStartOfInterval(loc);	
		char remoteBracket = loc.direction2 == BreakendDirection.Forward ? ']' : '[';
		String target = String.format("%s:%d", processContext.getDictionary().getSequence(loc.referenceIndex2).getSequenceName(), loc.start2);
		if (loc.direction == BreakendDirection.Forward) {
			alleles(referenceBase, String.format("%s%s%c%s%c", referenceBase, untemplatedSequence, remoteBracket, target, remoteBracket));
		} else {
			alleles(referenceBase, String.format("%c%s%c%s%s", remoteBracket, target, remoteBracket, untemplatedSequence, referenceBase));
		}
		return this;
	}
}

