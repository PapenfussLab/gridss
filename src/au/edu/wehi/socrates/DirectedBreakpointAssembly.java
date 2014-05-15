package au.edu.wehi.socrates;

import java.nio.charset.StandardCharsets;

import htsjdk.samtools.SAMRecord;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import au.edu.wehi.socrates.vcf.SvType;
import au.edu.wehi.socrates.vcf.VcfConstants;
import au.edu.wehi.socrates.vcf.VcfSvConstants;

public class DirectedBreakpointAssembly extends VariantContextDirectedBreakpoint {
	public final byte[] breakpointBaseQuality;
	public static DirectedBreakpointAssembly create(DirectedBreakpointAssembly variant, SAMRecord realignment) {
		if (realignment == null) return variant;
		VariantContextDirectedBreakpointBuilder builder = new VariantContextDirectedBreakpointBuilder(variant.processContext, variant)
			.realigned(variant, realignment);
		return new DirectedBreakpointAssembly(variant.processContext, builder.make());
	}
	public static DirectedBreakpointAssembly create(
			ProcessingContext processContext,
			String assemblerName,
			int referenceIndex,
			int position,
			BreakendDirection direction,
			byte[] breakpointSequence,
			byte[] breakpointBaseQuality,
			byte[] fullAssembly,
			byte[] fullAssemblyBaseQuality,
			int readCount,
			double breakpointQuality
			) {
		VariantContextBuilder builder = new VariantContextDirectedBreakpointBuilder(processContext)
			.location(new BreakendSummary(referenceIndex, direction, position, position, breakpointQuality), new String(breakpointSequence, StandardCharsets.US_ASCII))
			.id(String.format("%s-%s:%d-%s",
					assemblerName,
					processContext.getDictionary().getSequence(referenceIndex).getSequenceName(),
					position,
					direction == BreakendDirection.Forward ? "f" : "b"))
			.attribute(VcfConstants.ASSEMBLY_PROGRAM, assemblerName)
			.attribute(VcfConstants.ASSEMBLY_CONSENSUS, new String(fullAssembly, StandardCharsets.US_ASCII))
			.attribute(VcfConstants.ASSEMBLY_CONSENSUS_READ_COUNT, readCount);
		return new DirectedBreakpointAssembly(processContext, builder.make(), breakpointBaseQuality);
	}
	public static DirectedBreakpointAssembly create(
			ProcessingContext processContext,
			String assemblerName,
			int referenceIndex,
			int position,
			BreakendDirection direction,
			byte[] breakpointSequence,
			byte[] fullAssembly,
			Integer readCount,
			double breakpointQuality
			) {
		return create(processContext, assemblerName, referenceIndex, position, direction, breakpointSequence, null, fullAssembly, null, readCount, breakpointQuality);
	}
	protected DirectedBreakpointAssembly(ProcessingContext processContext, VariantContext variant) {
		this(processContext, variant, null);
	}
	protected DirectedBreakpointAssembly(ProcessingContext processContext, VariantContext variant, byte[] breakpointSequence) {
		super(processContext, variant);
		this.breakpointBaseQuality = breakpointSequence;
	}
	@Override
	public byte[] getBreakpointQuality() {
		return breakpointBaseQuality;
	}
	@Override
	public boolean isValid() {
		return super.isValid() && getAssemblerProgram() != null;
	};
	public String getAssemblerProgram() { return getAttributeAsString(VcfConstants.ASSEMBLY_PROGRAM, null); }
	public String getAssemblyConsensus() { return getAttributeAsString(VcfConstants.ASSEMBLY_CONSENSUS, ""); }
	public String getAnchorConsensus() {
		String consensus = getAssemblyConsensus();
		if (getBreakendSummary().direction == BreakendDirection.Forward) {
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
