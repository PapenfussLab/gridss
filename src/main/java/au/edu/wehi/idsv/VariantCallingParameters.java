package au.edu.wehi.idsv;

import htsjdk.variant.variantcontext.VariantContextBuilder;

import java.util.List;

import au.edu.wehi.idsv.vcf.VcfFilter;

import com.google.common.collect.Lists;

public class VariantCallingParameters {
	/**
	 * Minimum score for variant to be called
	 */
	public double minScore = 25;
	/**
	 * Maximum somatic p-value for flagging a variant as somatic
	 */
	public double somaticPvalueThreshold = 0.001;
	/**
	 * Call breakends only on assembled contigs
	 */
	public boolean callOnlyAssemblies = false;
	/**
	 * Minimum indel size
	 */
	public int minIndelSize = 16;
	/**
	 * Margin of error around breakends
	 * This margin is used to mitigate alignment errors around breakend coordinates
	 */
	public int breakendMargin = 10;
	public boolean writeFilteredCalls = Defaults.WRITE_FILTERED_CALLS;
	public BreakendSummary withMargin(ProcessingContext context, BreakendSummary bp) {
		return bp.expandBounds(breakendMargin);
	}
	public BreakendSummary withoutMargin(BreakendSummary bp) {
		return bp.compressBounds(breakendMargin);
	}
	public List<VcfFilter> breakpointFilters(BreakpointSummary bp) {
		List<VcfFilter> list = Lists.newArrayList();
		if (bp.couldBeDeletionOfSize(1, minIndelSize - 1)) {
			// likely to be an artifact
			// due to noise/poor alignment (eg bowtie2 2.1.0 would misalign reference reads)
			// and a nearby (real) indel
			// causing real indel mates to be assembled with noise read
			list.add(VcfFilter.SMALL_INDEL);
		}
		if (bp.couldBeReferenceAllele()) {
			list.add(VcfFilter.REFERENCE_ALLELE);
		}
		return list;
	}
	public VariantContextDirectedEvidence applyFilters(VariantContextDirectedEvidence call) {
		List<VcfFilter> filters = Lists.newArrayList();
		if (call instanceof VariantContextDirectedBreakpoint) {
			VariantContextDirectedBreakpoint vcdbp = (VariantContextDirectedBreakpoint)call;
			BreakpointSummary bp = vcdbp.getBreakendSummary();
			filters.addAll(breakpointFilters(bp));
			if (vcdbp.getBreakpointQual() < minScore || vcdbp.getBreakpointEvidenceCount(EvidenceSubset.ALL) == 0) {
				filters.add(VcfFilter.LOW_BREAKPOINT_SUPPORT);
			}
		}
		if (!filters.isEmpty()) {
			VariantContextBuilder builder = new VariantContextBuilder(call);
			for (VcfFilter f : filters) {
				builder.filter(f.filter());
			}
			call = (VariantContextDirectedEvidence)IdsvVariantContext.create(call.processContext, call.source, builder.make());
		}
		return call;
	}
}
