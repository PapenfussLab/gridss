package au.edu.wehi.socrates;

import java.nio.charset.StandardCharsets;

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
	public VariantContextDirectedBreakpointBuilder realigned(DirectedBreakpoint evidence, SAMRecord realignment) {		
		if (realignment.getReadUnmappedFlag()) {
			attribute(VcfConstants.REALIGNMENT_FAILURE, true);
		} else {
			BreakpointLocation localAnchor = evidence.getBreakpointLocation();
			BreakpointLocation remoteAnchor;
			int untemplatedSequenceLength;
			if ((localAnchor.direction == BreakpointDirection.Forward && realignment.getReadNegativeStrandFlag()) ||
					(localAnchor.direction == BreakpointDirection.Backward && !realignment.getReadNegativeStrandFlag())) {
				// negative strand template match means we flip the expected direction
				// since breakend sequence is on the +ve strand,
				// realigned forward breakends on +ve stand would indicate backward breakend
				// realigned backward breakends on +ve stand would indicate forward breakend
				remoteAnchor = new BreakpointLocation(realignment.getReferenceIndex(), BreakpointDirection.Forward,
						realignment.getAlignmentEnd(), realignment.getAlignmentEnd(), 0);
				untemplatedSequenceLength = SAMRecordUtil.getEndSoftClipLength(realignment);
			} else {
				// ACGT. => CGT breakpoint, -> ACGT.---.CGT
				remoteAnchor = new BreakpointLocation(realignment.getReferenceIndex(), BreakpointDirection.Backward,
						realignment.getAlignmentStart(), realignment.getAlignmentStart(), 0);
				untemplatedSequenceLength = SAMRecordUtil.getStartSoftClipLength(realignment);
			}			
			String untemplatedSequence = new String(evidence.getBreakpointSequence(), StandardCharsets.US_ASCII);
			if (localAnchor.direction == BreakpointDirection.Forward) {
				untemplatedSequence = untemplatedSequence.substring(0, untemplatedSequenceLength);
			} else {
				untemplatedSequence = untemplatedSequence.substring(untemplatedSequence.length() - untemplatedSequenceLength);
			}
			location(new BreakpointInterval(evidence.getBreakpointLocation(), remoteAnchor, localAnchor.qual), untemplatedSequence);
		}
		return this;
	}
	public VariantContextDirectedBreakpointBuilder breakpoint(DirectedEvidence evidence) {
		if (evidence instanceof DirectedBreakpoint) {
			DirectedBreakpoint bp = (DirectedBreakpoint)evidence;
			byte[] seq = bp.getBreakpointSequence();
			byte[] qual = bp.getBreakpointQuality();
			location(evidence.getBreakpointLocation(), new String(seq, StandardCharsets.US_ASCII));
		} else {
			location(evidence.getBreakpointLocation(), "");
		}
		return this;
	}
	/**
	 * Sets the VCF attributes for a breakend at the given location
	 * @param loc
	 * @return
	 */
	private String setVcfLocationAsStartOfInterval(BreakpointLocation loc) {
		this.chr(processContext.getReference().getSequenceDictionary().getSequence(loc.referenceIndex).getSequenceName());
		this.start(loc.start);
		this.stop(loc.start);
		if (loc.end != loc.start) {
			// TODO: correct setting of CI_POS for all directions
			this.attribute(VcfSvConstants.CONFIDENCE_INTERVAL_START_POSITION_KEY, 0);
			this.attribute(VcfSvConstants.CONFIDENCE_INTERVAL_END_POSITION_KEY, loc.end - loc.start);
		}
		String chr = processContext.getReference().getSequenceDictionary().getSequence(loc.referenceIndex).getSequenceName();
		return new String(processContext.getReference().getSubsequenceAt(chr, loc.start, loc.start).getBases(), StandardCharsets.US_ASCII);
	}
	public VariantContextDirectedBreakpointBuilder location(BreakpointLocation loc, String untemplatedSequence) {
		if (loc instanceof BreakpointInterval) {
			return location((BreakpointInterval)loc, "");
		}
		String referenceBase = setVcfLocationAsStartOfInterval(loc);
		String alt;
		if (loc.direction == BreakpointDirection.Forward) {
			alt = referenceBase + untemplatedSequence + '.';
		} else {
			alt = '.' + untemplatedSequence + referenceBase;
		}
		alleles(referenceBase, alt);
		log10PError(loc.qual / -10);
		return this;
	}
	public VariantContextDirectedBreakpointBuilder location(BreakpointInterval loc, String untemplatedSequence) {
		setVcfLocationAsStartOfInterval(loc);
		String referenceBase = setVcfLocationAsStartOfInterval(loc);	
		char remoteBracket = loc.direction2 == BreakpointDirection.Forward ? ']' : '[';
		String target = String.format("%s:%d", processContext.getDictionary().getSequence(loc.referenceIndex2).getSequenceName(), loc.start2);
		if (loc.direction == BreakpointDirection.Forward) {
			alleles(referenceBase, String.format("%s%s%c%s%c", referenceBase, untemplatedSequence, remoteBracket, target, remoteBracket));
		} else {
			alleles(referenceBase, String.format("%c%s%c%s%s", remoteBracket, target, remoteBracket, untemplatedSequence, referenceBase));
		}
		return this;
	}
}

