package au.edu.wehi.socrates;

import java.nio.charset.StandardCharsets;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;

import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;

import au.edu.wehi.socrates.vcf.SvType;
import au.edu.wehi.socrates.vcf.VcfConstants;
import au.edu.wehi.socrates.vcf.VcfSvConstants;

public class DirectedBreakpointAssembly extends VariantContextDirectedBreakpoint {
	public static final String SOURCE_NAME = "socrates";
	public static DirectedBreakpointAssembly create(DirectedBreakpointAssembly variant, SAMRecord realignment) {
		if (realignment == null) return variant;
		VariantContextBuilder builder = new VariantContextBuilder(variant);
		if (realignment.getReadUnmappedFlag()) {
			builder.attribute(VcfConstants.ASSEMBLY_REALIGNMENT_FAILURE, true);
			return new DirectedBreakpointAssembly(variant.processContext, variant);
		}
		BreakpointLocation remoteAnchor;
		int untemplatedSequenceLength;
		BreakpointDirection direction = variant.getBreakpointLocation().direction;
		char remoteBracket;
		if ((direction == BreakpointDirection.Forward && realignment.getReadNegativeStrandFlag()) ||
				(direction == BreakpointDirection.Backward && !realignment.getReadNegativeStrandFlag())) {
			// negative strand template match means we flip the expected direction
			// since breakend sequence is on the +ve strand,
			// realigned forward breakends on +ve stand would indicate backward breakend
			// realigned backward breakends on +ve stand would indicate forward breakend
			remoteAnchor = new BreakpointLocation(realignment.getReferenceIndex(), BreakpointDirection.Forward,
					realignment.getAlignmentEnd(), realignment.getAlignmentEnd(), 0);
			untemplatedSequenceLength = SAMRecordUtil.getEndSoftClipLength(realignment);
			remoteBracket = ']';
		} else {
			// ACGT. => CGT breakpoint, -> ACGT.---.CGT
			remoteAnchor = new BreakpointLocation(realignment.getReferenceIndex(), BreakpointDirection.Backward,
					realignment.getAlignmentStart(), realignment.getAlignmentStart(), 0);
			untemplatedSequenceLength = SAMRecordUtil.getStartSoftClipLength(realignment);
			remoteBracket = '[';
		}
		String untemplatedSequence = variant.getBreakpointSequenceString();
		String referenceBase = variant.getReference().getBaseString();
		String target = String.format("%s:%d", variant.processContext.getDictionary().getSequence(remoteAnchor.referenceIndex).getSequenceName(), remoteAnchor.start);
		if (direction == BreakpointDirection.Forward) {
			untemplatedSequence = untemplatedSequence.substring(0, untemplatedSequenceLength);
			builder.alleles(referenceBase, String.format("%s%s%c%s%c", referenceBase, untemplatedSequence, remoteBracket, target, remoteBracket));
		} else {
			untemplatedSequence = untemplatedSequence.substring(untemplatedSequence.length() - untemplatedSequenceLength);
			builder.alleles(referenceBase, String.format("%c%s%c%s%s", remoteBracket, target, remoteBracket, untemplatedSequence, referenceBase));
		}
		return new DirectedBreakpointAssembly(variant.processContext, builder.make());
	}
	public static DirectedBreakpointAssembly create(
			ProcessingContext processContext,
			String assemblerName,
			int referenceIndex,
			int position,
			BreakpointDirection direction,
			byte[] breakpointSequence,
			byte[] breakpointBaseQuality,
			byte[] fullAssembly,
			byte[] fullAssemblyBaseQuality,
			int readCount,
			double breakpointQuality
			) {
		return create(processContext, assemblerName, referenceIndex, position, direction, breakpointSequence, fullAssembly, readCount, breakpointQuality);
	}
	public static DirectedBreakpointAssembly create(
			ProcessingContext processContext,
			String assemblerName,
			int referenceIndex,
			int position,
			BreakpointDirection direction,
			byte[] breakpointSequence,
			byte[] fullAssembly,
			Integer readCount,
			double breakpointQuality
			) {
		String chr = processContext.getDictionary().getSequence(referenceIndex).getSequenceName();
		VariantContextBuilder builder = new VariantContextBuilder()
			.source(SOURCE_NAME)
			.id(String.format("%s-%s:%d-%s", assemblerName, chr, position, direction == BreakpointDirection.Forward ? "f" : "b"))
			.chr(chr)
			.start(position)
			.stop(position)
			.log10PError(breakpointQuality / -10)
			.attributes(null)
			.attribute(VcfSvConstants.SV_TYPE_KEY, SvType.BND.name())
			.attribute(VcfConstants.ASSEMBLY_PROGRAM, assemblerName)
			.attribute(VcfConstants.ASSEMBLY_CONSENSUS, new String(fullAssembly, StandardCharsets.US_ASCII))
			.attribute(VcfConstants.ASSEMBLY_CONSENSUS_READ_COUNT, readCount);
		String referenceBase = "N";
		if (processContext.getReference() != null) {
			referenceBase = new String(processContext.getReference().getSubsequenceAt(chr, position, position).getBases(), StandardCharsets.US_ASCII);
		}
		String alt;
		String breakStr = new String(breakpointSequence, StandardCharsets.US_ASCII);
		// TODO: use reference bases or anchor assembly?
		if (direction == BreakpointDirection.Forward) {
			alt = referenceBase + breakStr + '.';
		} else {
			alt = '.' + breakStr + referenceBase;
		}
		builder.alleles(referenceBase, alt);
		return new DirectedBreakpointAssembly(processContext, builder.make());
	}
	protected DirectedBreakpointAssembly(ProcessingContext processContext, VariantContext context) {
		super(processContext, context);
	}
	@Override
	public boolean isValid() {
		return super.isValid() && getAssemblerProgram() != null;
	};
	public String getAssemblerProgram() { return getAttributeAsString(VcfConstants.ASSEMBLY_PROGRAM, null); }
	public String getAssemblyConsensus() { return getAttributeAsString(VcfConstants.ASSEMBLY_CONSENSUS, ""); }
	public String getAnchorConsensus() {
		String consensus = getAssemblyConsensus();
		if (getBreakpointLocation().direction == BreakpointDirection.Forward) {
			return consensus.substring(0, consensus.length() - getBreakpointSequenceString().length());
		} else {
			return consensus.substring(getBreakpointSequenceString().length());
		}
	}
	/**
	 * Number of reads contributing to consensus
	 * @return
	 */
	public int getConsensusReadCount() {
		return getAttributeAsInt(VcfConstants.ASSEMBLY_CONSENSUS_READ_COUNT, 0);
	}
}
