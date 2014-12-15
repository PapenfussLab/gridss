package au.edu.wehi.idsv;

import htsjdk.variant.vcf.VCFConstants;

import java.nio.charset.StandardCharsets;
import java.util.Collections;
import java.util.List;

import org.apache.commons.lang3.StringUtils;

import au.edu.wehi.idsv.vcf.VcfAttributes;
import au.edu.wehi.idsv.vcf.VcfSvConstants;

import com.google.common.base.Function;
import com.google.common.collect.ComparisonChain;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Ordering;

public class StructuralVariationCallBuilder extends IdsvVariantContextBuilder {
	private final ProcessingContext processContext;
	private final VariantContextDirectedEvidence parent;
	private final List<SoftClipEvidence> scList = Lists.newArrayList();
	private final List<NonReferenceReadPair> rpList = Lists.newArrayList();
	private final List<AssemblyEvidence> assList = Lists.newArrayList();
	private BreakendSummary calledBreakend;
	public StructuralVariationCallBuilder(ProcessingContext processContext, VariantContextDirectedEvidence parent) {
		super(processContext, parent);
		this.processContext = processContext;
		this.parent = parent;
		calledBreakend = parent.getBreakendSummary();
	}
	public StructuralVariationCallBuilder addEvidence(DirectedEvidence evidence) {
		if (evidence == null) throw new NullPointerException();
		BreakendSummary bs = evidence.getBreakendSummary();
		if (evidence instanceof SoftClipEvidence) {
			scList.add((SoftClipEvidence)evidence);
			bs = processContext.getSoftClipParameters().withMargin(processContext, bs);
		} else if (evidence instanceof AssemblyEvidence) {
			assList.add((AssemblyEvidence)evidence);
		} else if (evidence instanceof NonReferenceReadPair) {
			rpList.add((NonReferenceReadPair)evidence);
		} else {
			throw new RuntimeException(String.format("Unknown evidence type %s", evidence.getClass()));
		}
		if (!parent.getBreakendSummary().overlaps(bs)) {
			throw new IllegalArgumentException(String.format("Sanity check failure: Evidence %s does not provide support for call at %s", evidence.getBreakendSummary(), parent.getBreakendSummary()));
		}
		return this;
	}
	public VariantContextDirectedEvidence make() {
		orderEvidence();
		setBreakpoint();
		setImprecise();
		aggregateAssemblyAttributes();
		aggregateReadPairAttributes();
		aggregateSoftClipAttributes();		
		setLlr();
		id(getID());
		VariantContextDirectedEvidence variant = (VariantContextDirectedEvidence)IdsvVariantContext.create(processContext, null, super.make());
		variant = calcSpv(variant);
		variant = processContext.getVariantCallingParameters().applyFilters(variant);
		return variant;
	}
	private void orderEvidence() {
		Collections.sort(assList, ByBestDesc);
		// Only count unique assemblies
		dedup(assList);
		Collections.sort(scList, ByBestDesc);
		Collections.sort(rpList, ByBestDesc);
	}
	private void dedup(List<? extends DirectedEvidence> list) {
		for (int i = list.size() - 1; i > 0; i--) {
			if (StringUtils.isNotEmpty(list.get(i).getEvidenceID()) && !list.get(i).getEvidenceID().equals(".")) {
				for (int j = i - 1; j >= 0; j--) {
					if (list.get(i).getEvidenceID().equals(list.get(j).getEvidenceID())) {
						list.remove(i);
						break;
					}
				}
			}
		}
	}
	private void setBreakpoint() {
		DirectedBreakpoint bestbp = null;
		if (assList.size() > 0 && assList.get(0) instanceof DirectedBreakpoint) {
			bestbp = (DirectedBreakpoint)assList.get(0);
		} else if (scList.size() > 0 && scList.get(0) instanceof DirectedBreakpoint) {
			bestbp = (DirectedBreakpoint)scList.get(0);
		}
		if (bestbp != null) {
			breakpoint(bestbp.getBreakendSummary(), bestbp.getUntemplatedSequence());
		}
	}
	private void setImprecise() {
		for (DirectedEvidence e : Iterables.concat(assList, scList, rpList)) {
			if (e instanceof DirectedBreakpoint && e.isBreakendExact()) {
				// clear flag
				rmAttribute(VcfSvConstants.IMPRECISE_KEY);
				return;
			}
		}
		attribute(VcfSvConstants.IMPRECISE_KEY, true);
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
		for (AssemblyEvidence e : assList) {
			assllr += Models.llr(e);
		}
		double rpllrn = AttributeConverter.asDoubleListOffset(parent.getAttribute(VcfAttributes.READPAIR_LOG_LIKELIHOOD_RATIO.attribute()), 0, 0d);
		double rpllrt = AttributeConverter.asDoubleListOffset(parent.getAttribute(VcfAttributes.READPAIR_LOG_LIKELIHOOD_RATIO.attribute()), 1, 0d);
		for (NonReferenceReadPair e : rpList) {
			if (e.getEvidenceSource().isTumour()) {
				rpllrt += Models.llr(e);
			} else {
				rpllrn += Models.llr(e);
			}
		}
		double scllrn = AttributeConverter.asDoubleListOffset(parent.getAttribute(VcfAttributes.SOFTCLIP_LOG_LIKELIHOOD_RATIO.attribute()), 0, 0d);
		double scllrt = AttributeConverter.asDoubleListOffset(parent.getAttribute(VcfAttributes.SOFTCLIP_LOG_LIKELIHOOD_RATIO.attribute()), 1, 0d);
		for (SoftClipEvidence e : scList) {
			if (e.getEvidenceSource().isTumour()) {
				scllrt += Models.llr(e);
			} else {
				scllrn += Models.llr(e);
			}
		}
		double beQual = 0;
		double bpQual = 0;
		for (DirectedEvidence e : Iterables.concat(assList, rpList, scList)) {
			if (e instanceof DirectedBreakpoint) {
				bpQual += Models.llr(e);
			} else {
				beQual += Models.llr(e);
			}
		}
		attribute(VcfAttributes.BREAKEND_QUAL.attribute(), beQual);
		attribute(VcfAttributes.BREAKPOINT_QUAL.attribute(), bpQual);
		attribute(VcfAttributes.CALLED_QUAL.attribute(), parent.getPhredScaledQual());
		attribute(VcfAttributes.ASSEMBLY_LOG_LIKELIHOOD_RATIO.attribute(), assllr);
		attribute(VcfAttributes.SOFTCLIP_LOG_LIKELIHOOD_RATIO.attribute(), ImmutableList.of(scllrn, scllrt ));
		attribute(VcfAttributes.READPAIR_LOG_LIKELIHOOD_RATIO.attribute(), ImmutableList.of(rpllrn, rpllrt ));
		double llr = assllr + scllrn + scllrt + rpllrn + rpllrt;
		attribute(VcfAttributes.LOG_LIKELIHOOD_RATIO.attribute(), llr);
		phredScore(llr);
	}
	private void aggregateAssemblyAttributes() {
		List<AssemblyEvidence> fullList = Lists.newArrayList(assList);
		List<AssemblyEvidence> supportList = Lists.newArrayList(fullList);
		attribute(VcfAttributes.ASSEMBLY_EVIDENCE_COUNT, fullList.size());
		attribute(VcfAttributes.ASSEMBLY_MAPPED, sumTN(supportList, new Function<AssemblyEvidence, Integer>() {
			@Override
			public Integer apply(AssemblyEvidence input) {
				return input instanceof DirectedBreakpoint ? 1 : 0;
			}
		}));
		attribute(VcfAttributes.ASSEMBLY_MAPQ_MAX, maxTN(supportList, new Function<AssemblyEvidence, Integer>() {
			@Override
			public Integer apply(AssemblyEvidence arg0) {
				return arg0.getLocalMapq();
			}}));
		attribute(VcfAttributes.ASSEMBLY_MAPQ_LOCAL_MAX, maxTN(supportList, new Function<AssemblyEvidence, Integer>() {
			@Override
			public Integer apply(AssemblyEvidence arg0) {
				return arg0.getLocalMapq();
			}}));
		attribute(VcfAttributes.ASSEMBLY_MAPQ_MAX, maxTN(supportList, new Function<AssemblyEvidence, Integer>() {
			@Override
			public Integer apply(AssemblyEvidence arg0) {
				if (arg0 instanceof DirectedBreakpoint) return ((DirectedBreakpoint)arg0).getRemoteMapq();
				return 0;
			}}));
		attribute(VcfAttributes.ASSEMBLY_MAPQ_TOTAL, sumTN(supportList, new Function<AssemblyEvidence, Integer>() {
			@Override
			public Integer apply(AssemblyEvidence arg0) {
				if (arg0 instanceof DirectedBreakpoint) return ((DirectedBreakpoint)arg0).getRemoteMapq();
				return 0;
			}}));
		attribute(VcfAttributes.ASSEMBLY_LENGTH_LOCAL_MAX, maxTN(supportList, new Function<AssemblyEvidence, Integer>() {
			@Override
			public Integer apply(AssemblyEvidence arg0) {
				return arg0.getAssemblyAnchorLength();
			}}));
		attribute(VcfAttributes.ASSEMBLY_LENGTH_REMOTE_MAX, maxTN(supportList, new Function<AssemblyEvidence, Integer>() {
			@Override
			public Integer apply(AssemblyEvidence arg0) {
				return arg0.getBreakendSequence().length;
			}}));
		int[] basecount = new int[2];
		int[] rpcount = new int[2];
		int[] sccount = new int[2];
		int[] sclen = new int[2];
		int[] rpmaxlen = new int[2];
		int[] scmaxlen = new int[2];
		List<String> consensus = Lists.newArrayList();
		for (AssemblyEvidence e : supportList) {
			basecount[0] += e.getAssemblyBaseCount(EvidenceSubset.NORMAL);
			basecount[1] += e.getAssemblyBaseCount(EvidenceSubset.TUMOUR);
			rpcount[0] += e.getAssemblySupportCountReadPair(EvidenceSubset.NORMAL);
			rpcount[1] += e.getAssemblySupportCountReadPair(EvidenceSubset.TUMOUR);
			sccount[0] += e.getAssemblySupportCountSoftClip(EvidenceSubset.NORMAL);
			sccount[1] += e.getAssemblySupportCountSoftClip(EvidenceSubset.TUMOUR);
			sclen[0] += e.getAssemblySoftClipLengthTotal(EvidenceSubset.NORMAL);
			sclen[1] += e.getAssemblySoftClipLengthTotal(EvidenceSubset.TUMOUR);
			rpmaxlen[0] = Math.max(rpmaxlen[0], e.getAssemblyReadPairLengthMax(EvidenceSubset.NORMAL));
			rpmaxlen[1] = Math.max(rpmaxlen[1], e.getAssemblyReadPairLengthMax(EvidenceSubset.TUMOUR));
			scmaxlen[0] = Math.max(scmaxlen[0], e.getAssemblySoftClipLengthMax(EvidenceSubset.NORMAL));
			scmaxlen[1] = Math.max(scmaxlen[1], e.getAssemblySoftClipLengthMax(EvidenceSubset.TUMOUR));
			consensus.add(new String(e.getAssemblySequence(), StandardCharsets.US_ASCII));
		}
		attribute(VcfAttributes.ASSEMBLY_BASE_COUNT, basecount);
		attribute(VcfAttributes.ASSEMBLY_READPAIR_COUNT, rpcount);
		attribute(VcfAttributes.ASSEMBLY_READPAIR_LENGTH_MAX, rpmaxlen);
		attribute(VcfAttributes.ASSEMBLY_SOFTCLIP_COUNT, sccount);
		attribute(VcfAttributes.ASSEMBLY_SOFTCLIP_CLIPLENGTH_TOTAL, sclen);
		attribute(VcfAttributes.ASSEMBLY_SOFTCLIP_CLIPLENGTH_MAX, scmaxlen);
		attribute(VcfAttributes.ASSEMBLY_CONSENSUS, consensus);
		
		if (calledBreakend instanceof BreakpointSummary) {
			int countrpn = 0;
			int countrpt = 0;
			int countscn = 0;
			int countsct = 0;
			for (AssemblyEvidence ae : supportList) {
				for (DirectedEvidence e : ae.getEvidence()) {
					if (e instanceof DirectedBreakpoint) {
						BreakpointSummary bp = ((DirectedBreakpoint)e).getBreakendSummary();
						if (!bp.overlaps(calledBreakend)) {
							if (e instanceof SoftClipEvidence) {
								if ((((SoftClipEvidence)e).getEvidenceSource()).isTumour()) {
									countsct++;
								} else {
									countscn++;
								}
							} else if (e instanceof NonReferenceReadPair) {
								if ((((NonReferenceReadPair)e).getEvidenceSource()).isTumour()) {
									countrpt++;
								} else {
									countrpn++;
								}
							}
						}
					}
				}
				attribute(VcfAttributes.ASSEMBLY_READPAIR_REMAPPED, new int[] { countrpn, countrpt });
				attribute(VcfAttributes.ASSEMBLY_SOFTCLIP_REMAPPED, new int[] { countscn, countsct });
			}
		}
	}
	private void aggregateReadPairAttributes() {
		attribute(VcfAttributes.READPAIR_EVIDENCE_COUNT, sumTN(rpList, new Function<NonReferenceReadPair, Integer>() {
			@Override
			public Integer apply(NonReferenceReadPair arg0) {
				return 1;
			}}));
		attribute(VcfAttributes.READPAIR_MAPPED_READPAIR, sumTN(rpList, new Function<NonReferenceReadPair, Integer>() {
			@Override
			public Integer apply(NonReferenceReadPair arg0) {
				return arg0 instanceof DirectedBreakpoint ? 1 : 0;
			}}));
		attribute(VcfAttributes.READPAIR_MAPQ_LOCAL_MAX, maxTN(rpList, new Function<NonReferenceReadPair, Integer>() {
			@Override
			public Integer apply(NonReferenceReadPair arg0) {
				return arg0.getLocalMapq();
			}}));
		attribute(VcfAttributes.READPAIR_MAPQ_LOCAL_TOTAL, sumTN(rpList, new Function<NonReferenceReadPair, Integer>() {
			@Override
			public Integer apply(NonReferenceReadPair arg0) {
				return arg0.getLocalMapq();
			}}));
		attribute(VcfAttributes.READPAIR_MAPQ_REMOTE_MAX, maxTN(rpList, new Function<NonReferenceReadPair, Integer>() {
			@Override
			public Integer apply(NonReferenceReadPair arg0) {
				if (arg0 instanceof DirectedBreakpoint) return ((DirectedBreakpoint)arg0).getRemoteMapq();
				return 0;
			}}));
		attribute(VcfAttributes.READPAIR_MAPQ_REMOTE_TOTAL, sumTN(rpList, new Function<NonReferenceReadPair, Integer>() {
			@Override
			public Integer apply(NonReferenceReadPair arg0) {
				if (arg0 instanceof DirectedBreakpoint) return ((DirectedBreakpoint)arg0).getRemoteMapq();
				return 0;
			}}));
	}
	private void aggregateSoftClipAttributes() {
		attribute(VcfAttributes.SOFTCLIP_EVIDENCE_COUNT, sumTN(scList, new Function<SoftClipEvidence, Integer>() {
			@Override
			public Integer apply(SoftClipEvidence arg0) {
				return 1;
			}}));
		attribute(VcfAttributes.SOFTCLIP_MAPPED, sumTN(scList, new Function<SoftClipEvidence, Integer>() {
			@Override
			public Integer apply(SoftClipEvidence arg0) {
				return arg0 instanceof RealignedSoftClipEvidence ? 1 : 0;
			}}));
		attribute(VcfAttributes.SOFTCLIP_MAPQ_REMOTE_TOTAL, sumTN(scList, new Function<SoftClipEvidence, Integer>() {
			@Override
			public Integer apply(SoftClipEvidence arg0) {
				if (arg0 instanceof DirectedBreakpoint) return ((DirectedBreakpoint)arg0).getRemoteMapq();
				return 0;
			}}));
		attribute(VcfAttributes.SOFTCLIP_MAPQ_REMOTE_MAX, maxTN(scList, new Function<SoftClipEvidence, Integer>() {
			@Override
			public Integer apply(SoftClipEvidence arg0) {
				if (arg0 instanceof DirectedBreakpoint) return ((DirectedBreakpoint)arg0).getRemoteMapq();
				return 0;
			}}));
		attribute(VcfAttributes.SOFTCLIP_LENGTH_REMOTE_TOTAL, sumTN(scList, new Function<SoftClipEvidence, Integer>() {
			@Override
			public Integer apply(SoftClipEvidence arg0) {
				return arg0.getSoftClipLength();
			}}));
		attribute(VcfAttributes.SOFTCLIP_LENGTH_REMOTE_MAX, maxTN(scList, new Function<SoftClipEvidence, Integer>() {
			@Override
			public Integer apply(SoftClipEvidence arg0) {
				return arg0.getSoftClipLength();
			}}));
	}
	private static <T extends DirectedEvidence> List<Integer> sumTN(Iterable<T> it, Function<T, Integer> f) {
		return sumOrMaxTN(it, f, false);
	}
	private static <T extends DirectedEvidence> List<Integer> maxTN(Iterable<T> it, Function<T, Integer> f) {
		return sumOrMaxTN(it, f, true);
	}
	private static <T extends DirectedEvidence> List<Integer> sumOrMaxTN(Iterable<T> it, Function<T, Integer> f, boolean max) {
		boolean collapse = false;
		List<Integer> list = Lists.newArrayList();
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
	private static Ordering<DirectedEvidence> ByBestDesc = new Ordering<DirectedEvidence>() {
		public int compare(DirectedEvidence o1, DirectedEvidence o2) {
			  return ComparisonChain.start()
			        .compareTrueFirst(o1 instanceof DirectedBreakpoint, o2 instanceof DirectedBreakpoint)
			        .compareTrueFirst(o1.isBreakendExact(), o2.isBreakendExact())
			        //.compareFalseFirst(o1 instanceof RemoteEvidence, o2 instanceof RemoteEvidence)
			        .compare(Models.llr(o2), Models.llr(o1)) // desc
			        .result();
		  }
	};
}
