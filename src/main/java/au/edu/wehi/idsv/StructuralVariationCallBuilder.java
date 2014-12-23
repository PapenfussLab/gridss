package au.edu.wehi.idsv;

import htsjdk.samtools.util.Log;
import htsjdk.variant.vcf.VCFConstants;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import au.edu.wehi.idsv.vcf.VcfAttributes;
import au.edu.wehi.idsv.vcf.VcfSvConstants;

import com.google.common.base.Predicate;
import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Ordering;

public class StructuralVariationCallBuilder extends IdsvVariantContextBuilder {
	private static final Log log = Log.getInstance(StructuralVariationCallBuilder.class);
	private final ProcessingContext processContext;
	private final VariantContextDirectedEvidence parent;
	private List<DirectedEvidence> list = Lists.newArrayList();
	public StructuralVariationCallBuilder(ProcessingContext processContext, VariantContextDirectedEvidence parent) {
		super(processContext, parent);
		this.processContext = processContext;
		this.parent = parent;
	}
	private BreakendSummary getBreakendWithMargin(DirectedEvidence evidence) {
		BreakendSummary bs = evidence.getBreakendSummary();
		if (evidence instanceof SoftClipEvidence) {
			bs = processContext.getSoftClipParameters().withMargin(processContext, bs);
		}
		return bs;
	}
	public StructuralVariationCallBuilder addEvidence(DirectedEvidence evidence) {
		if (evidence == null) throw new NullPointerException();
		list.add(evidence);
		BreakendSummary bs = getBreakendWithMargin(evidence);
		if (!parent.getBreakendSummary().overlaps(bs)) {
			throw new IllegalArgumentException(String.format("Sanity check failure: Evidence %s does not provide support for call at %s", evidence.getBreakendSummary(), parent.getBreakendSummary()));
		}
		return this;
	}
	public VariantContextDirectedEvidence make() {
		extractAssemblySupport();
		deduplicateEvidence();
		orderEvidence();
		setBreakpoint();
		setImprecise();
		attribute(VcfAttributes.CALLED_QUAL.attribute(), parent.getPhredScaledQual());
		setBreakpointAttributes();
		setBreakendAttributes();
		id(getID());
		VariantContextDirectedEvidence variant = (VariantContextDirectedEvidence)IdsvVariantContext.create(processContext, null, super.make());
		variant = calcSpv(variant);
		variant = processContext.getVariantCallingParameters().applyFilters(variant);
		return variant;
	}
	/**
	 * Includes the evidence supporting the assembly in the RP and SC totals 
	 */
	private void extractAssemblySupport() {
		for (AssemblyEvidence ae : Lists.newArrayList(Iterables.filter(list, AssemblyEvidence.class))) {
			list.addAll(ae.getEvidence());
		}
	}
	private void deduplicateEvidence() {
		Map<String, DirectedEvidence> unique = new HashMap<String, DirectedEvidence>();
		for (DirectedEvidence e : list) {
			unique.put(e.getEvidenceID(), e);
		}
		if (list.size() != unique.size()) {
			log.debug(String.format("Deduplicated %d records from %s", list.size() - unique.size(), getID()));
			list = new ArrayList<DirectedEvidence>(unique.values());
		}
	}
	private void orderEvidence() {
		Collections.sort(list, ByBestDesc);
	}
	private void setBreakpoint() {
		// Use breakpoint assembly, or breakpoint softclip to provide exact call bounds
		DirectedBreakpoint bestbp = (DirectedBreakpoint)Iterables.tryFind(list, new Predicate<DirectedEvidence>() {
			@Override
			public boolean apply(DirectedEvidence input) {
				return input.isBreakendExact() && input instanceof DirectedBreakpoint && input instanceof AssemblyEvidence;
			}
		}).or(Iterables.tryFind(list, new Predicate<DirectedEvidence>() {
			@Override
			public boolean apply(DirectedEvidence input) {
				return input.isBreakendExact() && input instanceof DirectedBreakpoint && input instanceof SoftClipEvidence;
			}
		})).orNull();
		if (bestbp != null) {
			breakpoint(bestbp.getBreakendSummary(), bestbp.getUntemplatedSequence());
		}
	}
	private void setImprecise() {
		if (Iterables.any(list, new Predicate<DirectedEvidence>() {
				@Override
				public boolean apply(DirectedEvidence input) {
					return input instanceof DirectedBreakpoint && input.isBreakendExact();
				}})) {
			rmAttribute(VcfSvConstants.IMPRECISE_KEY);
		} else {
			attribute(VcfSvConstants.IMPRECISE_KEY, true);
		}
	}
	private int offsetTN(DirectedEvidence e) {
		EvidenceSource source = e.getEvidenceSource();
		return (source instanceof SAMEvidenceSource && ((SAMEvidenceSource)source).isTumour()) ? 1 : 0;
	}
	private void setBreakendAttributes() {
		int asCount = 0;
		int[] scCount = new int[] {0, 0};
		int[] rpCount = new int[] {0, 0};
		float totalQual = 0;
		float asQual = 0;
		float[] scQual = new float[] {0, 0};
		float[] rpQual = new float[] {0, 0};
		for (DirectedEvidence e : list) {
			if (e instanceof DirectedBreakpoint) continue;
			float qual = e.getBreakendQual();
			int offset = offsetTN(e);
			if (e instanceof AssemblyEvidence) {
				if (e instanceof RemoteEvidence) {
					throw new RuntimeException("Sanity check failure: remote support for only breakend should not be possible");
					//rasCount++;
					//rasQual += qual;
				} else {
					asCount++;
					asQual += qual;
				}
			} else if (e instanceof SoftClipEvidence) {
				if (e instanceof RemoteEvidence) {
					throw new RuntimeException("Sanity check failure: remote support for only breakend should not be possible");
					//rscCount[offset]++;
					//rscQual[offset] += qual;
				} else {
					scCount[offset]++;
					scQual[offset] += qual;
				}
			} else if (e instanceof NonReferenceReadPair) {
				rpCount[offset]++;
				rpQual[offset] += qual;
			} else {
				throw new RuntimeException("Unknown evidence type " + e.getClass().toString());
			}
			totalQual += qual;
		}
		attribute(VcfAttributes.BREAKEND_ASSEMBLY_COUNT.attribute(), asCount); 
		attribute(VcfAttributes.BREAKEND_SOFTCLIP_COUNT.attribute(), scCount);
		attribute(VcfAttributes.BREAKEND_READPAIR_COUNT.attribute(), rpCount);
		attribute(VcfAttributes.BREAKEND_ASSEMBLY_QUAL.attribute(), asQual); 
		attribute(VcfAttributes.BREAKEND_SOFTCLIP_QUAL.attribute(), scQual);
		attribute(VcfAttributes.BREAKEND_READPAIR_QUAL.attribute(), rpQual);
		attribute(VcfAttributes.BREAKEND_QUAL.attribute(), totalQual);
	}
	private void setBreakpointAttributes() {
		int asCount = 0;
		int[] scCount = new int[] {0, 0};
		int[] rpCount = new int[] {0, 0};
		int rasCount = 0;
		int[] rscCount = new int[] {0, 0};
		float totalQual = 0;
		float asQual = 0;
		float[] scQual = new float[] {0, 0};
		float[] rpQual = new float[] {0, 0};
		float rasQual = 0;
		float[] rscQual = new float[] {0, 0};
		for (DirectedBreakpoint e : Iterables.filter(list, DirectedBreakpoint.class)) {
			float qual = e.getBreakpointQual();
			int offset = offsetTN(e);
			if (e instanceof AssemblyEvidence) {
				if (e instanceof RemoteEvidence) {
					rasCount++;
					rasQual += qual;
				} else {
					asCount++;
					asQual += qual;
				}
			} else if (e instanceof SoftClipEvidence) {
				if (e instanceof RemoteEvidence) {
					rscCount[offset]++;
					rscQual[offset] += qual;
				} else {
					scCount[offset]++;
					scQual[offset] += qual;
				}
			} else if (e instanceof NonReferenceReadPair) {
				rpCount[offset]++;
				rpQual[offset] += qual;
			} else {
				throw new RuntimeException("Unknown evidence type " + e.getClass().toString());
			}
			totalQual += qual;
		}
		attribute(VcfAttributes.BREAKPOINT_ASSEMBLY_COUNT.attribute(), asCount); 
		attribute(VcfAttributes.BREAKPOINT_SOFTCLIP_COUNT.attribute(), scCount);
		attribute(VcfAttributes.BREAKPOINT_READPAIR_COUNT.attribute(), rpCount);
		attribute(VcfAttributes.BREAKPOINT_ASSEMBLY_COUNT_REMOTE.attribute(), rasCount); 
		attribute(VcfAttributes.BREAKPOINT_SOFTCLIP_COUNT_REMOTE.attribute(), rscCount);		
		attribute(VcfAttributes.BREAKPOINT_ASSEMBLY_QUAL.attribute(), asQual); 
		attribute(VcfAttributes.BREAKPOINT_SOFTCLIP_QUAL.attribute(), scQual);
		attribute(VcfAttributes.BREAKPOINT_READPAIR_QUAL.attribute(), rpQual);
		attribute(VcfAttributes.BREAKPOINT_ASSEMBLY_QUAL_REMOTE.attribute(), rasQual); 
		attribute(VcfAttributes.BREAKPOINT_SOFTCLIP_QUAL_REMOTE.attribute(), rscQual);
		phredScore(totalQual);
	}
	private VariantContextDirectedEvidence calcSpv(VariantContextDirectedEvidence variant) {
		// TODO: somatic p-value should use evidence from both sides of the breakpoint
		float pvalue = Models.somaticPvalue(processContext, variant);
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(processContext, variant);
		builder.attribute(VcfAttributes.SOMATIC_P_VALUE, pvalue);
		if (pvalue < processContext.getVariantCallingParameters().somaticPvalueThreshold) {
			builder.attribute(VCFConstants.SOMATIC_KEY, true);
		} else {
			builder.rmAttribute(VCFConstants.SOMATIC_KEY);
		}
		return (VariantContextDirectedEvidence)builder.make();
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
			        .compare(o2 instanceof DirectedBreakpoint ? ((DirectedBreakpoint)o2).getBreakpointQual() : o2.getBreakendQual(),
			        		o1 instanceof DirectedBreakpoint ? ((DirectedBreakpoint)o1).getBreakpointQual() : o1.getBreakendQual()) // desc
			        .result();
		  }
	};
}
