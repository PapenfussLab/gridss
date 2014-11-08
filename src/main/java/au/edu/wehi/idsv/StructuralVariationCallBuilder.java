package au.edu.wehi.idsv;

import htsjdk.variant.vcf.VCFConstants;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.apache.commons.lang3.StringUtils;

import au.edu.wehi.idsv.vcf.VcfAttributes;
import au.edu.wehi.idsv.vcf.VcfAttributes.Subset;
import au.edu.wehi.idsv.vcf.VcfFilter;

import com.google.common.base.Function;
import com.google.common.base.Predicate;
import com.google.common.collect.ComparisonChain;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Ordering;

public class StructuralVariationCallBuilder extends IdsvVariantContextBuilder {
	private final ProcessingContext processContext;
	private final VariantContextDirectedEvidence parent;
	private final List<SoftClipEvidence> scList = new ArrayList<>();
	private final List<NonReferenceReadPair> rpList = new ArrayList<>();
	private final List<VariantContextDirectedEvidence> assList = new ArrayList<>();
	private BreakendSummary calledBreakend;
	public StructuralVariationCallBuilder(ProcessingContext processContext, VariantContextDirectedEvidence parent) {
		super(processContext, parent);
		this.processContext = processContext;
		this.parent = parent;
		calledBreakend = parent.getBreakendSummary();
	}
	public StructuralVariationCallBuilder addEvidence(DirectedEvidence evidence) {
		if (evidence == null) throw new NullPointerException();
		if (evidence instanceof SoftClipEvidence) {
			scList.add((SoftClipEvidence)evidence);
		} else if (evidence instanceof VariantContextDirectedEvidence) {
			assList.add((VariantContextDirectedEvidence)evidence);
		} else if (evidence instanceof NonReferenceReadPair) {
			rpList.add((NonReferenceReadPair)evidence);
		} else {
			throw new RuntimeException(String.format("Unknown evidence type %s", evidence.getClass()));
		}
		if (!parent.getBreakendSummary().overlaps(evidence.getBreakendSummary())) {
			throw new IllegalArgumentException(String.format("Sanity check failure: Evidence %s does not provide support for call at %s", evidence.getBreakendSummary(), parent.getBreakendSummary()));
		}
		return this;
	}
	public VariantContextDirectedEvidence make() {
		Collections.sort(assList, AssemblyByBestDesc);
		for (int i = assList.size() - 1; i > 0; i--) {
			// Only count unique assemblies
			if (StringUtils.isNotEmpty(assList.get(i).getID())
					&& !assList.get(i).getID().equals(".")
					&& assList.get(i).getID().equals(assList.get(i - 1).getID())) {
				assList.remove(i);
			}
		}
		Collections.sort(scList, SoftClipByBestDesc);
		VariantContextDirectedBreakpoint bestAssemblyBreakpoint = (VariantContextDirectedBreakpoint)Iterables.getFirst(Iterables.filter(assList, new Predicate<VariantContextDirectedEvidence>() {
			@Override
			public boolean apply(VariantContextDirectedEvidence input) {
				return input instanceof VariantContextDirectedBreakpoint
						// don't include remote breakpoints as their info (including coordinate) is of the remote position 
						&& !(input instanceof VariantContextDirectedBreakpointRemote);
			}
		}), null);
		if (bestAssemblyBreakpoint != null) {
			// override the breakpoint call to ensure untemplated sequence is present in the output
			calledBreakend = bestAssemblyBreakpoint.getBreakendSummary();
			breakpoint(bestAssemblyBreakpoint.getBreakendSummary(), bestAssemblyBreakpoint.getBreakpointSequenceString());
		} else if (scList.size() > 0 && scList.get(0) instanceof RealignedSoftClipEvidence) {
			RealignedSoftClipEvidence bp = (RealignedSoftClipEvidence)scList.get(0);
			calledBreakend = bp.getBreakendSummary();
			breakpoint(bp.getBreakendSummary(), bp.getUntemplatedSequence());
		}
		aggregateAssemblyAttributes();
		aggregateReadPairAttributes();
		aggregateSoftClipAttributes();		
		setLlr();
		id(getID());
		if (calledBreakend instanceof BreakpointSummary) {
			BreakpointSummary bp = (BreakpointSummary)calledBreakend;
			if (bp.referenceIndex == bp.referenceIndex2
					&& bp.direction != bp.direction2
					&& bp.end - bp.start == bp.end2 - bp.start2
					&& Math.abs(bp.start - bp.start2) < processContext.getVariantCallingParameters().minIndelSize
					) {
				// likely to be an artifact
				// due to noise/poor alignment (eg bowtie2 2.1.0 would misalign reference reads)
				// and a nearby (real) indel
				// causing real indel mates to be assembled with noise read
				filter(VcfFilter.SMALL_INDEL.filter());
			}
		}
		VariantContextDirectedEvidence variant = (VariantContextDirectedEvidence)IdsvVariantContext.create(processContext, null, super.make());
		variant = calcSpv(variant);
		return variant;
	}
	private VariantContextDirectedEvidence calcSpv(VariantContextDirectedEvidence variant) {
		// TODO: somatic p-value should use evidence from both sides of the breakpoint
		double pvalue = Models.somaticPvalue(processContext, variant);
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(processContext, variant);
		builder.attribute(VcfAttributes.SOMATIC_P_VALUE, pvalue);
		if (pvalue < processContext.getVariantCallingParameters().somaticPvalueThreshold) {
			builder.attribute(VCFConstants.SOMATIC_KEY, true);
		} else {
			builder.rmAttribute(VCFConstants.SOMATIC_KEY);
		}
		return (VariantContextDirectedEvidence)builder.make();
	}
	private void setLlr() {
		double assllr = parent.getAttributeAsDouble(VcfAttributes.ASSEMBLY_LOG_LIKELIHOOD_RATIO.attribute(), 0);
		for (VariantContextDirectedEvidence e : assList) {
			//assllr += PhredLogLikelihoodRatioModel.llr(e); // no need to recalculate as we already have the result stored
			assllr += e.getBreakendLogLikelihoodAssembly();
		}
		double rpllrn = parent.getAttributeAsDoubleListOffset(VcfAttributes.READPAIR_LOG_LIKELIHOOD_RATIO.attribute(), 0, 0d);
		double rpllrt = parent.getAttributeAsDoubleListOffset(VcfAttributes.READPAIR_LOG_LIKELIHOOD_RATIO.attribute(), 1, 0d);
		for (NonReferenceReadPair e : rpList) {
			if (e.getEvidenceSource().isTumour()) {
				rpllrt += Models.llr(e);
			} else {
				rpllrn += Models.llr(e);
			}
		}
		double scllrn = parent.getAttributeAsDoubleListOffset(VcfAttributes.SOFTCLIP_LOG_LIKELIHOOD_RATIO.attribute(), 0, 0d);
		double scllrt = parent.getAttributeAsDoubleListOffset(VcfAttributes.SOFTCLIP_LOG_LIKELIHOOD_RATIO.attribute(), 1, 0d);
		for (SoftClipEvidence e : scList) {
			if (e.getEvidenceSource().isTumour()) {
				scllrt += Models.llr(e);
			} else {
				scllrn += Models.llr(e);
			}
		}
		
		attribute(VcfAttributes.ASSEMBLY_LOG_LIKELIHOOD_RATIO.attribute(), assllr);
		attribute(VcfAttributes.SOFTCLIP_LOG_LIKELIHOOD_RATIO.attribute(), ImmutableList.of(scllrn, scllrt ));
		attribute(VcfAttributes.READPAIR_LOG_LIKELIHOOD_RATIO.attribute(), ImmutableList.of(rpllrn, rpllrt ));
		double llr = assllr + scllrn + scllrt + rpllrn + rpllrt;
		attribute(VcfAttributes.LOG_LIKELIHOOD_RATIO.attribute(), llr);
		phredScore(llr);
	}
	private List<Integer> setIntAttributeSumOrMax(final Iterable<VariantContextDirectedEvidence> list, final VcfAttributes attr, boolean max) {
		Iterable<List<Integer>> attrs = Iterables.transform(list, new Function<VariantContextDirectedEvidence, List<Integer>>() {
			@Override
			public List<Integer> apply(VariantContextDirectedEvidence arg0) {
				return arg0.getAttributeAsIntList(attr.attribute());
			}});
		List<Integer> result = new ArrayList<>();
		for (List<Integer> l : attrs) {
			while (l.size() > result.size()) result.add(0);
			for (int i = 0; i < l.size(); i++) {
				if (max) {
					result.set(i, Math.max(result.get(i), l.get(i)));
				} else {
					result.set(i, result.get(i) + l.get(i));
				}
			}
		}
		attribute(attr.attribute(), result);
		return result;
	}
	private List<String> setStringAttributeConcat(final Iterable<VariantContextDirectedEvidence> list, final VcfAttributes attr) {
		Iterable<List<String>> attrs = Iterables.transform(list, new Function<VariantContextDirectedEvidence, List<String>>() {
			@Override
			public List<String> apply(VariantContextDirectedEvidence arg0) {
				return arg0.getAttributeAsStringList(attr.attribute());
			}});
		List<String> result = Lists.newArrayList(Iterables.concat(attrs));
		attribute(attr.attribute(), result);
		return result;
	}
	private void aggregateAssemblyAttributes() {
		List<VariantContextDirectedEvidence> fullList = Lists.newArrayList(assList);
		List<VariantContextDirectedEvidence> list = Lists.newArrayList(fullList);
		if (processContext.getVariantCallingParameters().callOnlyAssemblies) {
			for (int i = list.size() - 1; i >= 0; i--) {
				BreakendSummary bs = list.get(i).getBreakendSummary();
				if (!(bs instanceof BreakpointSummary) || !((BreakpointSummary)bs).overlaps(calledBreakend)) {
					list.remove(i);
				}
			}
		}
		setIntAttributeSumOrMax(fullList, VcfAttributes.ASSEMBLY_EVIDENCE_COUNT, false);
		setIntAttributeSumOrMax(list, VcfAttributes.ASSEMBLY_MAPPED, false);
		setIntAttributeSumOrMax(list, VcfAttributes.ASSEMBLY_MAPQ_REMOTE_MAX, true);
		setIntAttributeSumOrMax(list, VcfAttributes.ASSEMBLY_MAPQ_REMOTE_TOTAL, false);
		setIntAttributeSumOrMax(list, VcfAttributes.ASSEMBLY_LENGTH_LOCAL_MAX, true);
		setIntAttributeSumOrMax(list, VcfAttributes.ASSEMBLY_LENGTH_REMOTE_MAX, true);
		setIntAttributeSumOrMax(list, VcfAttributes.ASSEMBLY_BASE_COUNT, false);
		setIntAttributeSumOrMax(list, VcfAttributes.ASSEMBLY_READPAIR_COUNT, false);
		setIntAttributeSumOrMax(list, VcfAttributes.ASSEMBLY_READPAIR_LENGTH_MAX, true);
		setIntAttributeSumOrMax(list, VcfAttributes.ASSEMBLY_SOFTCLIP_COUNT, false);
		setIntAttributeSumOrMax(list, VcfAttributes.ASSEMBLY_SOFTCLIP_CLIPLENGTH_TOTAL, false);
		setIntAttributeSumOrMax(list, VcfAttributes.ASSEMBLY_SOFTCLIP_CLIPLENGTH_MAX, true);
		setStringAttributeConcat(list, VcfAttributes.ASSEMBLY_CONSENSUS);
		setStringAttributeConcat(list, VcfAttributes.ASSEMBLY_PROGRAM);
		// setStringAttributeConcat(list, VcfAttributes.ASSEMBLY_BREAKEND_QUALS); // drop base quals
	}
	private void aggregateReadPairAttributes() {
		attribute(VcfAttributes.READPAIR_EVIDENCE_COUNT, sumOrMax(rpList, new Function<NonReferenceReadPair, Integer>() {
			@Override
			public Integer apply(NonReferenceReadPair arg0) {
				return 1;
			}}, false));
		attribute(VcfAttributes.READPAIR_MAPPED_READPAIR, sumOrMax(rpList, new Function<NonReferenceReadPair, Integer>() {
			@Override
			public Integer apply(NonReferenceReadPair arg0) {
				return arg0 instanceof DirectedBreakpoint ? 1 : 0;
			}}, false));
		attribute(VcfAttributes.READPAIR_MAPQ_LOCAL_MAX, sumOrMax(rpList, new Function<NonReferenceReadPair, Integer>() {
			@Override
			public Integer apply(NonReferenceReadPair arg0) {
				return arg0.getLocalMapq();
			}}, true));
		attribute(VcfAttributes.READPAIR_MAPQ_LOCAL_TOTAL, sumOrMax(rpList, new Function<NonReferenceReadPair, Integer>() {
			@Override
			public Integer apply(NonReferenceReadPair arg0) {
				return arg0.getLocalMapq();
			}}, false));
		attribute(VcfAttributes.READPAIR_MAPQ_REMOTE_MAX, sumOrMax(rpList, new Function<NonReferenceReadPair, Integer>() {
			@Override
			public Integer apply(NonReferenceReadPair arg0) {
				if (arg0 instanceof DirectedBreakpoint) return ((DirectedBreakpoint)arg0).getRemoteMapq();
				return 0;
			}}, true));
		attribute(VcfAttributes.READPAIR_MAPQ_REMOTE_TOTAL, sumOrMax(rpList, new Function<NonReferenceReadPair, Integer>() {
			@Override
			public Integer apply(NonReferenceReadPair arg0) {
				if (arg0 instanceof DirectedBreakpoint) return ((DirectedBreakpoint)arg0).getRemoteMapq();
				return 0;
			}}, false));
	}
	private void aggregateSoftClipAttributes() {
		attribute(VcfAttributes.SOFTCLIP_EVIDENCE_COUNT, sumOrMax(scList, new Function<SoftClipEvidence, Integer>() {
			@Override
			public Integer apply(SoftClipEvidence arg0) {
				return 1;
			}}, false));
		attribute(VcfAttributes.SOFTCLIP_MAPPED, sumOrMax(scList, new Function<SoftClipEvidence, Integer>() {
			@Override
			public Integer apply(SoftClipEvidence arg0) {
				return arg0 instanceof RealignedSoftClipEvidence ? 1 : 0;
			}}, false));
		attribute(VcfAttributes.SOFTCLIP_MAPQ_REMOTE_TOTAL, sumOrMax(scList, new Function<SoftClipEvidence, Integer>() {
			@Override
			public Integer apply(SoftClipEvidence arg0) {
				if (arg0 instanceof DirectedBreakpoint) return ((DirectedBreakpoint)arg0).getRemoteMapq();
				return 0;
			}}, false));
		attribute(VcfAttributes.SOFTCLIP_MAPQ_REMOTE_MAX, sumOrMax(scList, new Function<SoftClipEvidence, Integer>() {
			@Override
			public Integer apply(SoftClipEvidence arg0) {
				if (arg0 instanceof DirectedBreakpoint) return ((DirectedBreakpoint)arg0).getRemoteMapq();
				return 0;
			}}, true));
		attribute(VcfAttributes.SOFTCLIP_LENGTH_REMOTE_TOTAL, sumOrMax(scList, new Function<SoftClipEvidence, Integer>() {
			@Override
			public Integer apply(SoftClipEvidence arg0) {
				return arg0.getSoftClipLength();
			}}, false));
		attribute(VcfAttributes.SOFTCLIP_LENGTH_REMOTE_MAX, sumOrMax(scList, new Function<SoftClipEvidence, Integer>() {
			@Override
			public Integer apply(SoftClipEvidence arg0) {
				return arg0.getSoftClipLength();
			}}, true));
	}
	private static <T extends DirectedEvidence> List<Integer> sumOrMax(Iterable<T> it, Function<T, Integer> f, boolean max) {
		boolean collapse = false;
		List<Integer> list = new ArrayList<>();
		list.add(0);
		list.add(0);
		for (T t : it) {
			boolean isTumour = false;
			Integer v = f.apply(t);
			if (t.getEvidenceSource() instanceof SAMEvidenceSource) {
				isTumour = ((SAMEvidenceSource)t.getEvidenceSource()).isTumour();
			} else {
				collapse = true;
			}
			int offset = isTumour ? 1 : 0;
			if (max) {
				list.set(offset, Math.max(list.get(offset), v));
			} else {
				list.set(offset, list.get(offset) + v);
			}
			
		}
		if (collapse) {
			list.remove(1);
		}
		return list;
	}
	public StructuralVariationCallBuilder referenceReads(int normalCount, int tumourCount) {
		attribute(VcfAttributes.REFERENCE_COUNT_READ.attribute(), ImmutableList.of(normalCount, tumourCount ));
		return this;
	}
	public StructuralVariationCallBuilder referenceSpanningPairs(int normalCount, int tumourCount) {
		attribute(VcfAttributes.REFERENCE_COUNT_READPAIR.attribute(), ImmutableList.of(normalCount, tumourCount));
		return this;
	}
	private String getID() {
		BreakendSummary call = parent.getBreakendSummary();
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
	private static Ordering<SoftClipEvidence> SoftClipByBestDesc = new Ordering<SoftClipEvidence>() {
		public int compare(SoftClipEvidence o1, SoftClipEvidence o2) {
			  return ComparisonChain.start()
			        .compareTrueFirst(o1 instanceof RealignedSoftClipEvidence, o2 instanceof RealignedSoftClipEvidence)
			        .compare(Models.llr(o2), Models.llr(o2)) // desc
			        .result();
		  }
	};
	private static Ordering<VariantContextDirectedEvidence> AssemblyByBestDesc = new Ordering<VariantContextDirectedEvidence>() {
		public int compare(VariantContextDirectedEvidence o1, VariantContextDirectedEvidence o2) {
			  return ComparisonChain.start()
			        .compareTrueFirst(o1 instanceof VariantContextDirectedBreakpoint, o2 instanceof VariantContextDirectedBreakpoint)
			        .compareTrueFirst(o1.getAssemblySupportCountSoftClip(Subset.ALL) > 0, o2.getAssemblySupportCountSoftClip(Subset.ALL) > 0)
			        .compare(Models.llr(o2), Models.llr(o1)) // desc
			        .result();
		  }
	};
}
