package au.edu.wehi.idsv;

import au.edu.wehi.idsv.vcf.VcfFormatAttributes;
import gridss.analysis.InsertSizeDistribution;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class AlleleFractionAnnotator {
	private final ProcessingContext processContext;
	private final List<List<SAMEvidenceSource>> sourcesByCategory;
	private final int rpExclusionThreshold;
	public AlleleFractionAnnotator(ProcessingContext processContext, List<SAMEvidenceSource> sources) {
		this.processContext = processContext;
		this.sourcesByCategory = IntStream.of(processContext.getCategoryCount())
				.mapToObj(i -> sources.stream().filter(s -> s.getSourceCategory() == i).collect(Collectors.toList()))
				.collect(Collectors.toList());
		int maxFS = sources.stream().mapToInt(s -> s.getMaxConcordantFragmentSize()).max().orElse(500);
		int minFS = sources.stream().mapToInt(s -> s.getMinConcordantFragmentSize()).min().orElse(0);
		this.rpExclusionThreshold = maxFS + (maxFS - minFS);
	}

	public VariantContextDirectedEvidence annotate(VariantContextDirectedEvidence e) {
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(processContext, e);
		for (int category = 0; category < processContext.getCategoryCount(); category++) {
			builder.genotypeBuilder.get(category).attribute(VcfFormatAttributes.ALLELE_FRACTION.attribute(), calculateAF(e, category));
		}
		return (VariantContextDirectedEvidence)builder.make();
	}

	private double calculateAF(VariantContextDirectedEvidence e, int category) {
		int refCount = e.getReferenceReadCount(category);
		int varCount;
		boolean ignoreRP = false;
		if (e instanceof VariantContextDirectedBreakpoint) {
			VariantContextDirectedBreakpoint bp = (VariantContextDirectedBreakpoint) e;
			double expectedDP = 0;
			double expectedPortionAllocatedToReference = 0;
			BreakpointSummary bs = bp.getBreakendSummary();
			Integer eventSize = bs.getEventSize();
			if (eventSize != null && eventSize < rpExclusionThreshold) {
				// expectedDP = bp.getReferenceReadPairCount(category) * expectedPortionDiscordantPairUnderNullHypothesis(eventSize, category);
				// expectedPortionAllocatedToReference = expectedPortionOfVariantReadPairsAllocatedToReference(eventSize, category);
				ignoreRP = true;
			}
			// TODO: reduce varCount by expectedPortionDiscordantPairUnderNullHypothesis
			varCount = bp.getBreakpointEvidenceCount(category);
		} else {
			varCount = e.getBreakendEvidenceCount(category);
		}
		if (!ignoreRP) {
			refCount += e.getReferenceReadPairCount(category);
		}
		double totalSupport = refCount + varCount;
		// possible if all support is excluded RP support
		if (totalSupport == 0) return 0;
		return varCount  / totalSupport;
	}

	private double expectedPortionDiscordantPairUnderNullHypothesis(int deletionSize, int category) {
		return sourcesByCategory.get(category)
				.stream()
				.mapToDouble(s -> expectedPortionDiscordantPairUnderNullHypothesis(s, deletionSize))
				.max()
				.orElse(0);
	}

	/**
	 * We want to calculate the expected SV support under the null hypothesis.
	 * That is, what portion of the reference-support pairs we expect to indicate
	 * support for the (non-existent) variant.
	 */
	private double expectedPortionDiscordantPairUnderNullHypothesis(SAMEvidenceSource ses, int deletionSize) {
		InsertSizeDistribution isd = ses.getMetrics().getInsertSizeDistribution();
		if (isd == null) return 0;
		int maxFS = ses.getMaxConcordantFragmentSize();
		int minFS = ses.getMinConcordantFragmentSize();
		double concordantPortionOfDistribution = isd.cumulativeProbability(maxFS) - isd.cumulativeProbability(minFS);
		double expectedDiscordantPortion = isd.cumulativeProbability(maxFS + deletionSize);
		return expectedDiscordantPortion / concordantPortionOfDistribution;
	}
	private double expectedPortionOfVariantReadPairsAllocatedToReference(int deletionSize, int category) {
		return sourcesByCategory.get(category)
				.stream()
				.mapToDouble(s -> expectedPortionOfVariantReadPairsAllocatedToReference(s, deletionSize))
				.max()
				.orElse(0);
	}
	/**
	 * We want to know how badly we have messed up our RP estimate due to the distributions
	 * of the variant and reference overlapping. If our deletion is small then we'll have
	 * allocated a huge number of RPs that actually come from the variant to the reference.
	 */
	private double expectedPortionOfVariantReadPairsAllocatedToReference(SAMEvidenceSource ses, int deletionSize) {
		InsertSizeDistribution isd = ses.getMetrics().getInsertSizeDistribution();
		if (isd == null) return 0;
		int maxFS = ses.getMaxConcordantFragmentSize();
		int minFS = ses.getMinConcordantFragmentSize();
		double concordantPortionOfDistribution = isd.cumulativeProbability(maxFS) - isd.cumulativeProbability(minFS);
		double portionAllocatedToRef = isd.cumulativeProbability(maxFS - deletionSize) - isd.cumulativeProbability(minFS);
		return Math.max(portionAllocatedToRef / concordantPortionOfDistribution, 0);
	}
}
