package au.edu.wehi.socrates;

import java.util.List;

import com.google.common.collect.Lists;

public class StructuralVariationCallBuilder {
	private final ProcessingContext processContext;
	private final BreakpointLocation call;
	private final List<SoftClipEvidence> scList = Lists.newArrayList();
	private final List<NonReferenceReadPair> nrpList = Lists.newArrayList();
	private final List<DirectedBreakpointAssembly> assList = Lists.newArrayList();
	public StructuralVariationCallBuilder(ProcessingContext processContext, BreakpointLocation initialCall) {
		this.processContext = processContext;
		this.call  = initialCall;
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
		VariantContextDirectedBreakpointBuilder builder = createBuilder();
		builder.id(getID());
		return new StructuralVariationCall(processContext, builder.make());
	}
	private String getID() {
		
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
