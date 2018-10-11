package au.edu.wehi.idsv.configuration;

import java.util.Collection;
import java.util.HashSet;
import java.util.List;

import org.apache.commons.configuration.Configuration;

import com.google.common.collect.Lists;

import au.edu.wehi.idsv.BreakendSummary;
import au.edu.wehi.idsv.BreakpointSummary;
import au.edu.wehi.idsv.IdsvVariantContextBuilder;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;
import au.edu.wehi.idsv.VariantContextDirectedEvidence;
import au.edu.wehi.idsv.vcf.VcfFilter;

public class VariantCallingConfiguration {
	public static final String CONFIGURATION_PREFIX = "variantcalling";
	public VariantCallingConfiguration(Configuration config) {
		config = config.subset(CONFIGURATION_PREFIX);
		minReads = config.getDouble("minReads");
		minScore = config.getDouble("minScore");
		minSize = config.getInt("minSize");
		callUnassembledBreakpoints = config.getBoolean("callUnassembledBreakpoints");
		callUnassembledBreakends = config.getBoolean("callUnassembledBreakends");
		breakendMargin = config.getInt("breakendMargin");
		writeFiltered = config.getBoolean("writeFiltered");
		breakpointLowQuality = config.getDouble("breakpointLowQuality");
		breakendLowQuality = config.getDouble("breakendLowQuality");
		maxBreakendHomologyLength = config.getInt("maxBreakendHomologyLength");
		breakendHomologyAlignmentMargin = config.getInt("breakendHomologyAlignmentMargin");
		requireAssemblyCategorySupport = config.getBoolean("requireAssemblyCategorySupport");
		callBreakends = config.getBoolean("callBreakends");
		includeSupportingReadNames = config.getBoolean("includeSupportingReadNames");
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
	 * Minimum number of reads supporting variant either directly or indirectly through assembly
	 */
	public double minReads;
	/**
	 * Minimum score for variant to be called
	 */
	public double minScore;
	/**
	 * Minimum size of event for variant to be called
	 */
	public int minSize;
	/**
	 * Call only assembled breakpoint
	 */
	public boolean callUnassembledBreakpoints;
	/**
	 * Call only assembled breakends
	 */
	public boolean callUnassembledBreakends;
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
	public double breakpointLowQuality;
	public double breakendLowQuality;
	/**
	 * Maximum length of breakend homology to calculate
	 */
	public int maxBreakendHomologyLength;
	/**
	 * Number of reference bases to include in alignment
	 */
	public int breakendHomologyAlignmentMargin;
	/**
	 * Require that an anchored assembly is supported on both sides for each category
	 * Corrects for somatic calls with a flanking germline indel being called as somatic
	 * due to the non-zero germline support of the contig. 
	 */
	public boolean requireAssemblyCategorySupport;
	/**
	 * Include unpaired breakends in variant calls
	 */
	public boolean callBreakends;
	/**
	 * Include read names of all supporting reads
	 */
	public boolean includeSupportingReadNames;
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
	public List<VcfFilter> calculateCommonFilters(VariantContextDirectedEvidence call) {
		List<VcfFilter> filters = Lists.newArrayList();
		return filters;
	}
	public List<VcfFilter> calculateSingleBreakendFilters(VariantContextDirectedEvidence call) {
		List<VcfFilter> filters = Lists.newArrayList();
		if (!callUnassembledBreakends && call.getBreakendEvidenceCountAssembly() == 0) {
			filters.add(VcfFilter.NO_ASSEMBLY);
		}
		if (call.getBreakendQual() < minScore) {
			filters.add(VcfFilter.INSUFFICIENT_QUAL);
		}
		if (call.getBreakendSupportingFragmentCount() < minReads) {
			filters.add(VcfFilter.INSUFFICIENT_READS);
		}
		return filters;
	}
	public List<VcfFilter> calculateBreakpointFilters(VariantContextDirectedBreakpoint call) {
		List<VcfFilter> filters = Lists.newArrayList();
		BreakpointSummary bp = call.getBreakendSummary();
		if (call.getEventSize() != null && call.getEventSize() < minSize) {
			// over 90% of events are small. Since most SV analysis excludes such events
			// we allow the default output 
			filters.add(VcfFilter.SMALL_EVENT);
		}
		if (bp.couldBeReferenceAllele() && call.getUntemplatedSequence().length() == 0) {
			filters.add(VcfFilter.REFERENCE_ALLELE);
		}
		if (call.getBreakpointQual() < minScore) {
			filters.add(VcfFilter.INSUFFICIENT_QUAL);
		}
		if (call.getBreakpointSupportingFragmentCount() < minReads) {
			filters.add(VcfFilter.INSUFFICIENT_READS);
		}
		if (!callUnassembledBreakpoints && call.getBreakpointEvidenceCountAssembly() == 0) {
			filters.add(VcfFilter.NO_ASSEMBLY);
		}
		return filters;
	}
	public VariantContextDirectedEvidence applyConfidenceFilter(ProcessingContext processContext, final VariantContextDirectedEvidence variant) {
		Collection<String> filters = new HashSet<>(variant.getFilters());
		if (variant instanceof VariantContextDirectedBreakpoint) {
			VariantContextDirectedBreakpoint v = (VariantContextDirectedBreakpoint)variant;
			if (v.getBreakpointEvidenceCountLocalAssembly() == 0 && v.getBreakpointEvidenceCountRemoteAssembly() == 0) {
				filters.add(VcfFilter.NO_ASSEMBLY.filter());
			} else if (v.getBreakpointEvidenceCountLocalAssembly() == 0 || v.getBreakpointEvidenceCountRemoteAssembly() == 0) {
				filters.add(VcfFilter.SINGLE_ASSEMBLY.filter());
			} else if (v.getBreakpointEvidenceCountLocalAssembly() + v.getBreakpointEvidenceCountRemoteAssembly() > 0 &&
					v.getBreakpointEvidenceCountReadPair() + v.getBreakpointEvidenceCountSoftClip() == 0) {
				filters.add(VcfFilter.ASSEMBLY_ONLY.filter());
			}
			if (variant.getPhredScaledQual() < processContext.getVariantCallingParameters().breakpointLowQuality) {
				filters.add(VcfFilter.LOW_QUAL.filter());
			}
		} else {
			if (variant.getBreakendEvidenceCountAssembly() == 0) {
				filters.add(VcfFilter.NO_ASSEMBLY.filter());
			}
			if (variant.getPhredScaledQual() < processContext.getVariantCallingParameters().breakendLowQuality) {
				filters.add(VcfFilter.LOW_QUAL.filter());
			}
		}
		VariantContextDirectedEvidence filteredVariant = (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(processContext, variant).filters(filters.toArray(new String[0])).make();
		return filteredVariant;
	}
}
