package au.edu.wehi.idsv;

import htsjdk.variant.variantcontext.VariantContextBuilder;

import java.util.List;

import au.edu.wehi.idsv.vcf.VcfFilter;

import com.google.common.collect.Lists;

public class VariantCallingParameters {
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
	public int minIndelSize = 100; // DREAM challenge SV min size is 100bp //16;
	public List<VcfFilter> breakpointFilters(BreakpointSummary bp) {
		List<VcfFilter> list = Lists.newArrayList();
		if (bp.referenceIndex == bp.referenceIndex2
				&& bp.direction != bp.direction2
				&& bp.end - bp.start == bp.end2 - bp.start2
				&& Math.abs(bp.start - bp.start2) < minIndelSize
				) {
			// likely to be an artifact
			// due to noise/poor alignment (eg bowtie2 2.1.0 would misalign reference reads)
			// and a nearby (real) indel
			// causing real indel mates to be assembled with noise read
			list.add(VcfFilter.SMALL_INDEL);
		}
		return list;
	}
	public VariantContextDirectedEvidence applyFilters(VariantContextDirectedEvidence call) {
		List<VcfFilter> filters = Lists.newArrayList();
		if (call.getBreakendSummary() instanceof BreakpointSummary) {
			BreakpointSummary bp = (BreakpointSummary)call.getBreakendSummary();
			filters.addAll(breakpointFilters(bp));
		}
		if (call.getEvidenceCount(EvidenceSubset.ALL) == 0) {
			filters.add(VcfFilter.NO_SUPPORT);
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
