package au.edu.wehi.idsv.configuration;

import java.util.List;

import org.apache.commons.configuration.Configuration;

import com.google.common.collect.Lists;

import au.edu.wehi.idsv.BreakendSummary;
import au.edu.wehi.idsv.BreakpointSummary;
import au.edu.wehi.idsv.EvidenceSubset;
import au.edu.wehi.idsv.IdsvVariantContextBuilder;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;
import au.edu.wehi.idsv.VariantContextDirectedEvidence;
import au.edu.wehi.idsv.vcf.VcfFilter;

public class VariantCallingConfiguration {
	public static final String CONFIGURATION_PREFIX = "variantcalling";
	public VariantCallingConfiguration(Configuration config) {
		config = config.subset(CONFIGURATION_PREFIX);
		minScore = config.getDouble("minScore");
		callOnlyAssemblies = config.getBoolean("callOnlyAssemblies");
		minIndelSize = config.getInt("minIndelSize");
		breakendMargin = config.getInt("breakendMargin");
		writeFiltered = config.getBoolean("writeFiltered");
		switch (config.getString("format")) {
			case "vcf4.2":
				placeholderBreakend = false;
				break;
			case "vcf4.1":
				placeholderBreakend = true;
				break;
			default:
				throw new IllegalArgumentException(String.format("Unrecognised output format \"%d\"", config.getString("format")));
		}
	}
	/**
	 * Minimum score for variant to be called
	 */
	public double minScore;
	/**
	 * Call breakends only on assembled contigs
	 */
	public boolean callOnlyAssemblies;
	/**
	 * Minimum indel size
	 */
	public int minIndelSize;
	/**
	 * Number bases in which nearby evidence will be considered to support the same variant.
	 * This margin is used to mitigate soft clip alignment errors and microhomologies around breakend coordinates
	 */
	public int breakendMargin;
	public boolean writeFiltered;
	public boolean placeholderBreakend;
	public BreakendSummary withMargin(BreakendSummary bp) {
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
	public VariantContextDirectedEvidence applyConfidenceFilter(ProcessingContext processContext, final VariantContextDirectedEvidence variant) {
		VariantContextDirectedEvidence filteredVariant = variant;
		if (variant instanceof VariantContextDirectedBreakpoint) {
			VariantContextDirectedBreakpoint v = (VariantContextDirectedBreakpoint)filteredVariant;
			if (v.getBreakpointEvidenceCountLocalAssembly() == 0 && v.getBreakpointEvidenceCountRemoteAssembly() == 0) { 
				filteredVariant = (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(processContext, filteredVariant).filter(VcfFilter.NO_ASSEMBLY.filter()).make();
			} else if (v.getBreakpointEvidenceCountLocalAssembly() == 0 || v.getBreakpointEvidenceCountRemoteAssembly() == 0) {
				filteredVariant = (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(processContext, filteredVariant).filter(VcfFilter.SINGLE_ASSEMBLY.filter()).make();
			}
		} else {
			if (filteredVariant.getBreakendEvidenceCountAssembly() == 0) {
				filteredVariant = (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(processContext, filteredVariant).filter(VcfFilter.NO_ASSEMBLY.filter()).make();
			}
		}
		return filteredVariant;
	}
}
