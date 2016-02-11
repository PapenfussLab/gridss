package au.edu.wehi.idsv.configuration;

import java.util.List;

import org.apache.commons.configuration.Configuration;

import au.edu.wehi.idsv.BreakendSummary;
import au.edu.wehi.idsv.BreakpointSummary;
import au.edu.wehi.idsv.IdsvVariantContextBuilder;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;
import au.edu.wehi.idsv.VariantContextDirectedEvidence;
import au.edu.wehi.idsv.vcf.VcfFilter;

import com.google.common.collect.Lists;

public class VariantCallingConfiguration {
	public static final String CONFIGURATION_PREFIX = "variantcalling";
	public VariantCallingConfiguration(Configuration config) {
		config = config.subset(CONFIGURATION_PREFIX);
		minScore = config.getDouble("minScore");
		callOnlyAssemblies = config.getBoolean("callOnlyAssemblies");
		breakendMargin = config.getInt("breakendMargin");
		writeFiltered = config.getBoolean("writeFiltered");
		lowQuality = config.getDouble("lowQuality");
//		switch (config.getString("format")) {
//			case "vcf4.2":
//				placeholderBreakend = false;
//				break;
//			case "vcf4.1":
//				placeholderBreakend = true;
//				break;
//			default:
//				throw new IllegalArgumentException(String.format("Unrecognised output format \"%s\"", config.getString("format")));
//		}
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
	 * Number bases in which nearby evidence will be considered to support the same variant.
	 * This margin is used to mitigate soft clip alignment errors and microhomologies around breakend coordinates
	 */
	public int breakendMargin;
	/**
	 * Breakpoint event size (in multiples of breakendMargin) at which the full margin is applied
	 */
	private int fullMarginMultiple = 2;
	public boolean writeFiltered;
	//public boolean placeholderBreakend;
	public double lowQuality;
	public BreakendSummary withMargin(BreakendSummary bp) {
		if (bp == null) return null;
		return bp.expandBounds(marginFor(bp));
	}
	private int marginFor(BreakendSummary be) {
		if (be instanceof BreakpointSummary) {
			BreakpointSummary bp = (BreakpointSummary) be;
			if (bp.referenceIndex == bp.referenceIndex2) {
				int minsize = getMinSize(bp);
				if (minsize < breakendMargin * fullMarginMultiple) {
					return minsize / fullMarginMultiple;
				}
			}
		}
		return breakendMargin;
	}
	public int getMinSize(BreakpointSummary bp) {
		int minsize  = Math.min(Math.abs(bp.start - bp.end2), Math.abs(bp.start2 - bp.end));
		if (bp.localBreakend().overlaps(bp.remoteBreakend())) {
			minsize = 0;
		}
		return minsize;
	}
	public List<VcfFilter> calculateBreakendFilters(VariantContextDirectedEvidence call) {
		List<VcfFilter> filters = Lists.newArrayList();
		return filters;
	}
	public List<VcfFilter> calculateBreakpointFilters(VariantContextDirectedBreakpoint call) {
		List<VcfFilter> filters = Lists.newArrayList();
		BreakpointSummary bp = call.getBreakendSummary();
		
		//if (minIndelSize > 0 && bp.couldBeDeletionOfSize(1, minIndelSize - 1)) {
			// likely to be an artifact
			// due to noise/poor alignment (eg bowtie2 2.1.0 would misalign reference reads)
			// and a nearby (real) indel
			// causing real indel mates to be assembled with noise read
			//filters.add(VcfFilter.SMALL_INDEL);
		//}
		if (bp.couldBeReferenceAllele() && call.getUntemplatedSequence().length() == 0) {
			filters.add(VcfFilter.REFERENCE_ALLELE);
		}
		if (call.getBreakpointQual() < minScore || call.getBreakpointEvidenceCount() == 0) {
			filters.add(VcfFilter.LOW_BREAKPOINT_SUPPORT);
		}
		if (call.getBreakpointEvidenceCountAssembly() == 0 && call.getBreakpointEvidenceCountReadPair() + call.getBreakpointEvidenceCountSoftClip() == 1) {
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
			} else if (v.getBreakpointEvidenceCountLocalAssembly() + v.getBreakpointEvidenceCountRemoteAssembly() > 0 &&
					v.getBreakpointEvidenceCountReadPair() + v.getBreakpointEvidenceCountSoftClip() == 0) {
				filteredVariant = (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(processContext, filteredVariant).filter(VcfFilter.ASSEMBLY_ONLY.filter()).make();
			}
		} else {
			if (filteredVariant.getBreakendEvidenceCountAssembly() == 0) {
				filteredVariant = (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(processContext, filteredVariant).filter(VcfFilter.NO_ASSEMBLY.filter()).make();
			}
		}
		if (filteredVariant.getPhredScaledQual() < processContext.getVariantCallingParameters().lowQuality) {
			filteredVariant = (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(processContext, filteredVariant).filter(VcfFilter.LOW_QUAL.filter()).make();
		}
		return filteredVariant;
	}
}
