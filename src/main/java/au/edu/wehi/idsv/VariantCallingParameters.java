package au.edu.wehi.idsv;

import htsjdk.variant.variantcontext.VariantContextBuilder;

import java.util.List;

import au.edu.wehi.idsv.vcf.VcfFilter;

import com.google.common.collect.Lists;

public class VariantCallingParameters {
	/**
	 * Minimum score for variant to be called
	 */
	public double minScore = 38;
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
	/**
	 * Maximum coverage before evidence is filtered
	 */
	public int maxCoverage = 100000;
	public boolean writeFilteredCalls = Defaults.WRITE_FILTERED_CALLS;
	public BreakendSummary withMargin(ProcessingContext context, BreakendSummary bp) {
		if (bp == null) return null;
		return bp.expandBounds(breakendMargin);
	}
	public BreakendSummary withoutMargin(BreakendSummary bp) {
		return bp.compressBounds(breakendMargin);
	}
	public List<VcfFilter> calculateBreakendFilters(VariantContextDirectedEvidence call) {
		List<VcfFilter> filters = Lists.newArrayList();
		return filters;
	}
	public List<VcfFilter> calculateBreakpointFilters(VariantContextDirectedBreakpoint call) {
		List<VcfFilter> filters = Lists.newArrayList();
		BreakpointSummary bp = call.getBreakendSummary();
		if (minIndelSize > 0 && bp.couldBeDeletionOfSize(1, minIndelSize - 1)) {
			// likely to be an artifact
			// due to noise/poor alignment (eg bowtie2 2.1.0 would misalign reference reads)
			// and a nearby (real) indel
			// causing real indel mates to be assembled with noise read
			filters.add(VcfFilter.SMALL_INDEL);
		}
		if (bp.couldBeReferenceAllele() && call.getUntemplatedSequence().length() == 0) {
			filters.add(VcfFilter.REFERENCE_ALLELE);
		}
		if (call.getBreakpointQual() < minScore || call.getBreakpointEvidenceCount(EvidenceSubset.ALL) == 0) {
			filters.add(VcfFilter.LOW_BREAKPOINT_SUPPORT);
		}
		if (call.getBreakpointEvidenceCountAssembly() == 0 && call.getBreakpointEvidenceCountReadPair(EvidenceSubset.ALL) + call.getBreakpointEvidenceCountSoftClip(EvidenceSubset.ALL) == 1) {
			filters.add(VcfFilter.SINGLE_SUPPORT);
		}
		return filters;
	}
	public VariantContextDirectedEvidence applyConfidenceFilter(VariantContextDirectedEvidence variant) {
		VariantContextDirectedEvidence filteredVariant = variant;
		if (variant instanceof VariantContextDirectedBreakpoint) {
			VariantContextDirectedBreakpoint v = (VariantContextDirectedBreakpoint)filteredVariant;
			if (v.getBreakpointEvidenceCountLocalAssembly() == 0 && v.getBreakpointEvidenceCountRemoteAssembly() == 0) { 
				filteredVariant = (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(filteredVariant.processContext, filteredVariant).filter(VcfFilter.NO_ASSEMBLY.filter()).make();
			} else if (v.getBreakpointEvidenceCountLocalAssembly() == 0 || v.getBreakpointEvidenceCountRemoteAssembly() == 0) {
				variant = (VariantContextDirectedEvidence)new VariantContextBuilder(variant).filter(VcfFilter.SINGLE_ASSEMBLY.filter()).make();
			}
		} else {
			if (filteredVariant.getBreakendEvidenceCountAssembly() == 0) {
				filteredVariant = (VariantContextDirectedEvidence)new VariantContextBuilder(filteredVariant).filter(VcfFilter.NO_ASSEMBLY.filter()).make();
			}
		}
		return filteredVariant;
	}
}
