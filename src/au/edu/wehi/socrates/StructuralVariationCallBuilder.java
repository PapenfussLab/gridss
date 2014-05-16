package au.edu.wehi.socrates;

import htsjdk.variant.variantcontext.VariantContextBuilder;

import java.util.List;

import au.edu.wehi.socrates.vcf.VcfConstants;
import au.edu.wehi.socrates.vcf.VcfSvConstants;

import com.google.common.collect.Lists;

public class StructuralVariationCallBuilder {
	private final ProcessingContext processContext;
	private final BreakendSummary call;
	private final List<SoftClipEvidence> scList = Lists.newArrayList();
	private final List<NonReferenceReadPair> nrpList = Lists.newArrayList();
	private final List<VariantContextDirectedBreakpoint> assList = Lists.newArrayList();
	private int referenceReadCount = -1;
	private int referenceSpanningPairCount = -1;
	public StructuralVariationCallBuilder(ProcessingContext processContext, BreakendSummary call) {
		this.processContext = processContext;
		this.call  = call;
	}
	public StructuralVariationCallBuilder evidence(DirectedEvidence evidence) {
		if (evidence == null) throw new NullPointerException();
		if (evidence instanceof SoftClipEvidence) {
			scList.add((SoftClipEvidence)evidence);
		} else if (evidence instanceof VariantContextDirectedBreakpoint) {
			assList.add((VariantContextDirectedBreakpoint)evidence);
		} else if (evidence instanceof NonReferenceReadPair) {
			nrpList.add((NonReferenceReadPair)evidence);
		} else {
			throw new RuntimeException(String.format("Unknown evidence type %s", evidence.getClass()));
		}
		return this;
	}
	public VariantContextDirectedBreakpoint make() {
		VariantContextBuilder builder = createBuilder()
			.id(getID())
			.log10PError(call.evidence.getScore() / -10); // qual = phred for now
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
		sb.append(call.direction == BreakendDirection.Forward ? 'f' : 'b');
		if (call instanceof BreakpointSummary) {
			BreakpointSummary loc = (BreakpointSummary)call;
			sb.append(processContext.getDictionary().getSequence(loc.referenceIndex2).getSequenceName());
			sb.append(':');
			sb.append(loc.start2);
			if (loc.end2 != loc.start2) {
				sb.append('-');
				sb.append(loc.end2);
			}
			sb.append(loc.direction2 == BreakendDirection.Forward ? 'f' : 'b');
		}
		return sb.toString();
	}
	private VariantContextDirectedBreakpointBuilder createBuilder() {
		VariantContextDirectedBreakpointBuilder builder;
		VariantContextDirectedBreakpoint ass = bestAssembly();
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
			if (sce.getBreakendSummary() instanceof BreakpointSummary
					&& !(best.getBreakendSummary() instanceof BreakpointSummary)) {
				best = sce;
			} else if (sce.getBreakendSummary().getClass() == best.getBreakendSummary().getClass()) {
				if (sce.getBreakendSummary() instanceof BreakpointSummary) {
					// both BreakpointSummary
					if (sce.getSoftClipRealignmentSAMRecord().getMappingQuality() > best.getSoftClipRealignmentSAMRecord().getMappingQuality()) {
						// we have a better soft-clip mapping
						best = sce;
					}
				} else {
					// both BreakendSummary
					if (sce.getSoftClipLength() > best.getSoftClipLength()) {
						best = sce;
					}
				}
			}
		}
		return best;
	}
	private VariantContextDirectedBreakpoint bestAssembly() {
		if (assList.size() == 0) return null;
		VariantContextDirectedBreakpoint best = assList.get(0);
		for (VariantContextDirectedBreakpoint a : assList) {
			if (a.getPhredScaledQual() > best.getPhredScaledQual()) {
				best = a;
			}
		}
		return best;
	}
}
