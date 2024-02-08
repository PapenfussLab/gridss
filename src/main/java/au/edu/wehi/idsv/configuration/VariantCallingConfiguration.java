package au.edu.wehi.idsv.configuration;

import au.edu.wehi.idsv.*;
import au.edu.wehi.idsv.vcf.VcfFilter;
import au.edu.wehi.idsv.vcf.VcfSvConstants;
import com.google.common.collect.Lists;
import org.apache.commons.configuration.Configuration;

import java.util.Collection;
import java.util.HashSet;
import java.util.List;

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
		callBreakends = config.getBoolean("callBreakends");
		includeSupportingReadNames = config.getBoolean("includeSupportingReadNames");
		breakendMaxAssemblySupportBias = config.getDouble("breakendMaxAssemblySupportBias");
		callFullyAnchoredAssemblyVariants = config.getBoolean("callFullyAnchoredAssemblyVariants");
		ignoreMissingAssemblyFile = config.getBoolean("ignoreMissingAssemblyFile");
		minimumImpreciseDeletion = config.getInt("minimumImpreciseDeletion");
		requireReadPair = config.getBoolean("requireReadPair");
		requireSplitRead = config.getBoolean("requireSplitRead");
		includeSingleAssemblyFilter = config.getBoolean("includeSingleAssemblyFilter");
		requiredReadAndAssemblyBreakpointOverlap = config.getInt("requiredReadAndAssemblyBreakpointOverlap");
	}
	/**
	 * Ignore missing assembly file
	 */
	public boolean ignoreMissingAssemblyFile;
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
	 * Include unpaired breakends in variant calls
	 */
	public boolean callBreakends;
	/**
	 * Include read names of all supporting reads
	 */
	public boolean includeSupportingReadNames;
	/**
	 * Single breakend maximum support bias.
	 * Bias of 0 indicates the same number of directly supporting reads and read supporting via assembly (typically these are the same reads, although this is not strictly necessary)
	 * Bias of -1 indicates no assembly support
	 * Bias of 1 indicates no direct read support
	 */
	private double breakendMaxAssemblySupportBias;
	/**
	 * Determine whether include assembly evidence when the assembled variant is composed entirely of
	 * anchored assembly bases.
	 *
	 * We exclude these as the positional de Bruijn graph based GRIDSS assembler does not restrict assembly to
	 * reads contributing to this breakend. This means we may have reads in our anchored assembly that occur on
	 * a different haplotype to the reads assembled into the breakend potion of the contig
	 */
	public boolean callFullyAnchoredAssemblyVariants;
	/**
	 * Minimum size to call an imprecise deletion
	 * This filters out false positives caused by the reference-supporting reads at the
	 * edge of the library fragment size distribution.
	 */
	public int minimumImpreciseDeletion;
	public boolean requireReadPair;
	public boolean requireSplitRead;
	/**
	 * Include single_assembly filtering in output vcf
	 */
	public boolean includeSingleAssemblyFilter;
	/**
	 * Minimum number of bases that consider overlap between read and the breakpoint that the assembly generates
	 */
	public int requiredReadAndAssemblyBreakpointOverlap;
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
		if (bp.couldBeReferenceAllele(call.isBreakendExact()) && call.getUntemplatedSequence().length() == 0) {
			// unhandled edge case: precise only on one side
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
			if (!callUnassembledBreakpoints && (v.getBreakpointEvidenceCountLocalAssembly() == 0 && v.getBreakpointEvidenceCountRemoteAssembly() == 0) ){
				filters.add(VcfFilter.NO_ASSEMBLY.filter());
			} else if (includeSingleAssemblyFilter && (v.getBreakpointEvidenceCountLocalAssembly() == 0 || v.getBreakpointEvidenceCountRemoteAssembly() == 0)) {
				filters.add(VcfFilter.SINGLE_ASSEMBLY.filter());
			} else if (v.getBreakpointEvidenceCountLocalAssembly() + v.getBreakpointEvidenceCountRemoteAssembly() + v.getBreakpointEvidenceCountCompoundAssembly() > 0 &&
					v.getBreakpointEvidenceCountReadPair() + v.getBreakpointEvidenceCountSoftClip() + v.getBreakpointEvidenceCountIndel() == 0) {
				filters.add(VcfFilter.ASSEMBLY_ONLY.filter());
			}
			if (variant.getPhredScaledQual() < processContext.getVariantCallingParameters().breakpointLowQuality) {
				filters.add(VcfFilter.LOW_QUAL.filter());
			}
		} else {
			if (!callUnassembledBreakends && (variant.getBreakendEvidenceCountAssembly() == 0) ){
				filters.add(VcfFilter.NO_ASSEMBLY.filter());
			}
			if (variant.getPhredScaledQual() < processContext.getVariantCallingParameters().breakendLowQuality) {
				filters.add(VcfFilter.LOW_QUAL.filter());
			}
			// more filters
			int bassr = variant.getBreakendEvidenceCountAssemblySoftClip();
			int basrp = variant.getBreakendEvidenceCountAssemblyReadPair();
			int bsc = variant.getBreakendEvidenceCountSoftClip();
			int bum = variant.getBreakendEvidenceCountReadPair();
			int asmSupport = bassr + basrp;
			int directSupport = bsc + bum;
			double bebias = (asmSupport - directSupport) / (double)(asmSupport + directSupport);
			if (Math.abs(bebias) > breakendMaxAssemblySupportBias) {
				filters.add(VcfFilter.ASSEMBLY_BIAS.filter());
			}
			if ( (bsc == 0) && (requireSplitRead) ) {
				filters.add(VcfFilter.NO_SR.filter());
			}
			if ((bum == 0) && (requireReadPair)) {
				filters.add(VcfFilter.NO_RP.filter());
			}
		}
		VariantContextDirectedEvidence filteredVariant = (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(processContext, variant).filters(filters.toArray(new String[0])).make();
		return filteredVariant;
	}

	/**
	 * Determines whether to filter this putative variant before
	 * evidence allocation or
	 * Note that the variant call is minimally annotated at this point.
	 * Only the position, quality score (CQ), IMPRECISE flag is populated at this point.
	 * @param minimalVariant minimally annotated variant
	 * @return true if the variant should be ignored. False otherwise
	 */
    public boolean isHardFilteredBeforeAnnotation(final VariantContextDirectedEvidence minimalVariant) {
		// If we're under min score with all possible evidence allocated, we're definitely going to fail
		// when we restrict evidence to single breakpoint support
		if (minimalVariant.getBreakendQual() < minScore) return true;
		if (minimalVariant.getBreakendSummary() instanceof BreakpointSummary) {
			boolean isImprecise = minimalVariant.hasAttribute(VcfSvConstants.IMPRECISE_KEY);
			BreakpointSummary bp = ((BreakpointSummary)minimalVariant.getBreakendSummary());
			// Drop FP calls from fragments just outside the expected library fragment size bounds
			if (isImprecise && bp.couldBeDeletionOfSize(0, minimumImpreciseDeletion)) return true;
		}
		return false;
    }
}
