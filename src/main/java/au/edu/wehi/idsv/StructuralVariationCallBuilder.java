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
	public StructuralVariationCallBuilder addEvidence(DirectedEvidence evidence) {
		if (evidence == null) throw new NullPointerException();
		if (!isSupportingEvidence(evidence)) {
			throw new IllegalArgumentException(String.format("Sanity check failure: Evidence %s does not provide support for call at %s", evidence.getBreakendSummary(), parent.getBreakendSummary()));
		}
		list.add(evidence);
		return this;
	}
	private BreakendSummary getBreakendSummaryWithMargin(DirectedEvidence evidence) {
		BreakendSummary bs = evidence.getBreakendSummary();
		bs = processContext.getVariantCallingParameters().withMargin(processContext, bs);
		return bs;
	}
	private boolean isSupportingEvidence(DirectedEvidence evidence) {
		BreakendSummary bs = getBreakendSummaryWithMargin(evidence);
		return parent.getBreakendSummary().overlaps(bs);
	}
	public VariantContextDirectedEvidence make() {
		//extractAssemblySupport(); // don't extract - as extraction inflates evidence beyond original call
		deduplicateEvidence();
		orderEvidence();
		setBreakpoint();
		setImprecise();
		attribute(VcfAttributes.CALLED_QUAL.attribute(), parent.getPhredScaledQual());
		setBreakpointAttributes();
		setBreakendAttributes();
		// id(parent.getID()); // can't change from parent ID as the id is already referenced in the MATEID of the other breakend  
		VariantContextDirectedEvidence variant = (VariantContextDirectedEvidence)IdsvVariantContext.create(processContext, null, super.make());
		variant = calcSpv(variant);
		variant = processContext.getVariantCallingParameters().applyFilters(variant);
		//assert(sanitycheck(variant));
		return variant;
	}
	private boolean sanitycheck(VariantContextDirectedEvidence annotated) {
		double qual = annotated.getPhredScaledQual();
		double origqual = parent.getPhredScaledQual();
		assert(qual <= origqual + 0.01);
		return true;
	}
	/**
	 * Includes the evidence supporting the assembly in the RP and SC totals
	 * 
	 * ISSUE: would like to include alternatively mapped evidence as support but
	 * a) assemblies sometimes include read not even supporting the breakend
	 * b) extracting causes qual to be greater than called qual which causes problems with evidence assignment
	 */
	@SuppressWarnings("unused")
	private void extractAssemblySupport() {
		for (AssemblyEvidence ae : Lists.newArrayList(Iterables.filter(list, AssemblyEvidence.class))) {
			for (DirectedEvidence e : ae.getEvidence()) {
				if (ae instanceof DirectedBreakpoint) {
					list.add(e);
				} else {
					if (e instanceof DirectedBreakpoint && !isSupportingEvidence(e)) {
						// don't include breakpoint evidence from a breakend assembly
					} else {
						list.add(e);
					}
				}
			}
		}
	}
	private void deduplicateEvidence() {
		Map<String, DirectedEvidence> unique = new HashMap<String, DirectedEvidence>();
		for (DirectedEvidence e : list) {
			unique.put(e.getEvidenceID(), e);
		}
		if (list.size() != unique.size()) {
			log.debug(String.format("Deduplicated %d records from %s", list.size() - unique.size(), parent.getID()));
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
		String homo = "";
		if (bestbp != null) {
			breakpoint(bestbp.getBreakendSummary(), bestbp.getUntemplatedSequence());
			homo = bestbp.getHomologySequence();
		}
		if (homo.length() > 0) {
			attribute(VcfSvConstants.HOMOLOGY_SEQUENCE_KEY, homo);
			attribute(VcfSvConstants.HOMOLOGY_LENGTH_KEY, homo.length());
		} else {
			rmAttribute(VcfSvConstants.HOMOLOGY_SEQUENCE_KEY);
			rmAttribute(VcfSvConstants.HOMOLOGY_LENGTH_KEY);
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
		if (source instanceof SAMEvidenceSource) return ((SAMEvidenceSource)source).getSourceCategory();
		return 0;
	}
	private void setBreakendAttributes() {
		int asCount = 0;
		int[] scCount = new int[processContext.getCategoryCount()];
		int[] rpCount = new int[processContext.getCategoryCount()];
		float totalQual = 0;
		float asQual = 0;
		float[] scQual = new float[processContext.getCategoryCount()];
		float[] rpQual = new float[processContext.getCategoryCount()];
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
		int[] scCount = new int[processContext.getCategoryCount()];
		int[] rpCount = new int[processContext.getCategoryCount()];
		int rasCount = 0;
		int[] rscCount = new int[processContext.getCategoryCount()];
		float totalQual = 0;
		float asQual = 0;
		float[] scQual = new float[processContext.getCategoryCount()];
		float[] rpQual = new float[processContext.getCategoryCount()];
		float rasQual = 0;
		float[] rscQual = new float[processContext.getCategoryCount()];
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
	private static Ordering<DirectedEvidence> ByBestDesc = new Ordering<DirectedEvidence>() {
		public int compare(DirectedEvidence o1, DirectedEvidence o2) {
			ComparisonChain chain = ComparisonChain.start()
			        .compareTrueFirst(o1 instanceof DirectedBreakpoint, o2 instanceof DirectedBreakpoint)
			        .compareTrueFirst(o1.isBreakendExact(), o2.isBreakendExact())
			        .compare(o2 instanceof DirectedBreakpoint ? ((DirectedBreakpoint)o2).getBreakpointQual() : o2.getBreakendQual(),
			        		o1 instanceof DirectedBreakpoint ? ((DirectedBreakpoint)o1).getBreakpointQual() : o1.getBreakendQual()) // desc
			        .compareTrueFirst(o1 instanceof AssemblyEvidence, o2 instanceof AssemblyEvidence);
			if (o1 instanceof DirectedBreakpoint && o2 instanceof DirectedBreakpoint) {
				// TODO: favour small indels
				chain = chain.compare(((DirectedBreakpoint)o1).getBreakendSummary(), ((DirectedBreakpoint)o2).getBreakendSummary(), BreakpointSummary.ByLowHigh);
			}
			return chain.result();
		  }
	};
}
