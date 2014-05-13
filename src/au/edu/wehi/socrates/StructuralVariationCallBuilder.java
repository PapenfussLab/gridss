package au.edu.wehi.socrates;

import htsjdk.variant.variantcontext.VariantContextBuilder;

import java.util.List;

import au.edu.wehi.socrates.vcf.VcfConstants;
import au.edu.wehi.socrates.vcf.VcfSvConstants;

import com.google.common.collect.Lists;

public class StructuralVariationCallBuilder {
	private final ProcessingContext processContext;
	private final BreakpointLocation call;
	private final List<SoftClipEvidence> scList = Lists.newArrayList();
	private final List<NonReferenceReadPair> nrpList = Lists.newArrayList();
	private final List<DirectedBreakpointAssembly> assList = Lists.newArrayList();
	private int referenceReadCount = -1;
	private int referenceSpanningPairCount = -1;
	public StructuralVariationCallBuilder(ProcessingContext processContext, BreakpointLocation call) {
		this.processContext = processContext;
		this.call  = call;
	}
	public StructuralVariationCallBuilder evidence(DirectedEvidence evidence) {
		if (evidence == null) throw new NullPointerException();
		if (evidence instanceof SoftClipEvidence) {
			scList.add((SoftClipEvidence)evidence);
		} else if (evidence instanceof DirectedBreakpointAssembly) {
			assList.add((DirectedBreakpointAssembly)evidence);
		} else if (evidence instanceof NonReferenceReadPair) {
			nrpList.add((NonReferenceReadPair)evidence);
		} else {
			throw new RuntimeException(String.format("Unknown evidence type %s", evidence.getClass()));
		}
		return this;
	}
	public StructuralVariationCall make() {
		VariantContextBuilder builder = createBuilder()
			.id(getID())
			.log10PError(call.qual / -10); // qual = phred for now
		int oea = 0, dp = 0;
		for (NonReferenceReadPair nrp : nrpList) {
			if (nrp.getRemoteReferenceIndex() == -1) {
				oea++;
			} else {
				dp++;
			}
		}
		builder
			.attribute(VcfConstants.SOFT_CLIP_READ_COUNT, scList.size())
			.attribute(VcfConstants.DISCORDANT_READ_PAIR_COUNT, dp)
			.attribute(VcfConstants.UNMAPPED_MATE_READ_COUNT, oea);
		
		if (referenceReadCount >= 0) {
			builder.attribute(VcfConstants.REFERENCE_READ_COUNT, referenceReadCount);
		}
		if (referenceSpanningPairCount >= 0) {
			builder.attribute(VcfConstants.REFERENCE_SPANNING_READ_PAIR_COUNT, referenceSpanningPairCount);
		}
		return new StructuralVariationCall(processContext, builder.make());
	}
	public StructuralVariationCallBuilder referenceReads(int count) {
		referenceReadCount = count;
		return this;
	}
	public StructuralVariationCallBuilder referenceSpanningPairs(int count) {
		referenceSpanningPairCount = count;
		return this;
	}
	private String getID() {
		StringBuilder sb = new StringBuilder("call");
		sb.append(processContext.getDictionary().getSequence(call.referenceIndex).getSequenceName());
		sb.append(':');
		sb.append(call.start);
		if (call.end != call.start) {
			sb.append('-');
			sb.append(call.end);
		}
		sb.append(call.direction == BreakpointDirection.Forward ? 'f' : 'b');
		if (call instanceof BreakpointInterval) {
			BreakpointInterval loc = (BreakpointInterval)call;
			sb.append(processContext.getDictionary().getSequence(loc.referenceIndex2).getSequenceName());
			sb.append(':');
			sb.append(loc.start2);
			if (loc.end2 != loc.start2) {
				sb.append('-');
				sb.append(loc.end2);
			}
			sb.append(loc.direction2 == BreakpointDirection.Forward ? 'f' : 'b');
		}
		return sb.toString();
	}
	private VariantContextDirectedBreakpointBuilder createBuilder() {
		VariantContextDirectedBreakpointBuilder builder;
		DirectedBreakpointAssembly ass = bestAssembly();
		SoftClipEvidence sce = bestSoftclip();
		if (ass != null) {
			builder = new VariantContextDirectedBreakpointBuilder(processContext, ass);
		} else if (sce != null) {
			builder = new VariantContextDirectedBreakpointBuilder(processContext)
				.breakpoint(sce);
			if (sce.getSoftClipRealignmentSAMRecord() != null) {
				builder.realigned(sce, sce.getSoftClipRealignmentSAMRecord());
			}
		} else {
			builder = new VariantContextDirectedBreakpointBuilder(processContext)
				.location(call, "");
		}
		return builder;
	}
	private SoftClipEvidence bestSoftclip() {
		if (scList.size() == 0) return null;
		SoftClipEvidence best = scList.get(0);
		for (SoftClipEvidence sce : scList) {
			if (sce.getBreakpointLocation() instanceof BreakpointInterval
					&& !(best.getBreakpointLocation() instanceof BreakpointInterval)) {
				best = sce;
			} else if (sce.getBreakpointLocation().getClass() == best.getBreakpointLocation().getClass()) {
				if (sce.getBreakpointLocation() instanceof BreakpointInterval) {
					// both BreakpointInterval
					if (sce.getSoftClipRealignmentSAMRecord().getMappingQuality() > best.getSoftClipRealignmentSAMRecord().getMappingQuality()) {
						// we have a better soft-clip mapping
						best = sce;
					}
				} else {
					// both BreakpointLocation
					if (sce.getSoftClipLength() > best.getSoftClipLength()) {
						best = sce;
					}
				}
			}
		}
		return best;
	}
	private DirectedBreakpointAssembly bestAssembly() {
		if (assList.size() == 0) return null;
		DirectedBreakpointAssembly best = assList.get(0);
		for (DirectedBreakpointAssembly a : assList) {
			if (a.getPhredScaledQual() > best.getPhredScaledQual()) {
				best = a;
			}
		}
		return best;
	}
}
