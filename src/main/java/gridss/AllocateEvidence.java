package gridss;

import java.util.Iterator;
import java.util.concurrent.ExecutorService;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import com.google.common.collect.Iterators;

import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.DirectedEvidenceOrder;
import au.edu.wehi.idsv.GreedyVariantAllocationCache;
import au.edu.wehi.idsv.SAMEvidenceSource;
import au.edu.wehi.idsv.SequentialEvidenceAllocator;
import au.edu.wehi.idsv.SequentialEvidenceAllocator.VariantEvidenceSupport;
import au.edu.wehi.idsv.StructuralVariationCallBuilder;
import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;
import au.edu.wehi.idsv.configuration.VariantCallingConfiguration;
import au.edu.wehi.idsv.util.AsyncBufferedIterator;
import au.edu.wehi.idsv.util.AutoClosingIterator;
import au.edu.wehi.idsv.validation.OrderAssertingIterator;
import au.edu.wehi.idsv.validation.PairedEvidenceTracker;
import gridss.cmdline.VcfTransformCommandLineProgram;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;

@CommandLineProgramProperties(
        summary = "Evidence reallocation is required to ensure that any given read/read pair/assembly, "
        		+ " supports only a consistent set breakpoints. \n"
        		+ "For discordant read pairs, this means that only a single breakpoint can be supported. \n"
        		+ "For split reads, only a single set of partial read mappings can be supported, and only one breakpoint per split can be supported.\n"
        		+ "For indels, only a single read mapping can be supported, and only one breakpoint per indel can be supported.\n"
        		+ "Note: the EID attribute must be populated with the relevant GRIDSS EvidenceID's using AllocateEvidence",  
        oneLineSummary = "Uniquely allocates evidence supporting multiple mutually exclusive breakpoints.",
        programGroup=picard.cmdline.programgroups.Metrics.class
)
public class AllocateEvidence extends VcfTransformCommandLineProgram {
	private static final Log log = Log.getInstance(AllocateEvidence.class);
	@Argument(doc="Evidence allocation strategy used to uniquely assign evidence.")
	public EvidenceAllocationStrategy ALLOCATION_STRATEGY = EvidenceAllocationStrategy.GREEDY;
	public static enum EvidenceAllocationStrategy {
		GREEDY,
	}
	private GreedyVariantAllocationCache cache;
	public CloseableIterator<DirectedEvidence> getEvidenceIterator() {
		CloseableIterator<DirectedEvidence> evidenceIt;
		boolean assemblyOnly = getContext().getVariantCallingParameters().callOnlyAssemblies;
		if (assemblyOnly) {
			evidenceIt = SAMEvidenceSource.mergedIterator(ImmutableList.of(getAssemblySource()), true);
		} else {
			evidenceIt = SAMEvidenceSource.mergedIterator(ImmutableList.<SAMEvidenceSource>builder().addAll(getSamEvidenceSources()).add(getAssemblySource()).build(), true);
		}
		if (Defaults.SANITY_CHECK_ITERATORS) {
			evidenceIt = new AutoClosingIterator<>(
					new PairedEvidenceTracker<>("Evidence",
							new OrderAssertingIterator<>(evidenceIt, DirectedEvidenceOrder.ByNatural)),
					evidenceIt);
		}
		return evidenceIt;
	}
	@Override
	public CloseableIterator<VariantContextDirectedBreakpoint> iterator(CloseableIterator<VariantContextDirectedBreakpoint> calls, ExecutorService threadpool) {
		boolean multimapping = Iterables.any(getSamEvidenceSources(), ses -> ses.getMetrics().getIdsvMetrics().SECONDARY_NOT_SPLIT > 0);
		if (multimapping) {
			if (!getContext().getConfig().multimapping) {
				log.info("Not performing unique variant allocation of multimapping reads as the configuration setting multimapping is set to false.");
			} else if (!getContext().getConfig().multimappingUniqueVariantAllocation) {
				log.info("Not performing unique variant allocation of multimapping reads as the configuration setting multimappingUniqueVariantAllocation is set to false.");
			} else {
				log.info("Multimapping mode invoked due to existence of at least one BAM file with a non-split secondary alignment.");
				populateCache();
			}
		}
		log.info("Allocating evidence"); 
		CloseableIterator<DirectedEvidence> evidence = new AsyncBufferedIterator<>(getEvidenceIterator(), "mergedEvidence-allocation");
		Iterator<VariantEvidenceSupport> annotator = new SequentialEvidenceAllocator(getContext(), calls, evidence, SAMEvidenceSource.maximumWindowSize(getContext(), getSamEvidenceSources(), getAssemblySource()), true);
		Iterator<VariantContextDirectedBreakpoint> it = Iterators.transform(annotator, bp -> annotate(bp));
		it = Iterators.filter(it, v -> v != null);
		return new AutoClosingIterator<>(it, calls, evidence, cache);
	}
	public void populateCache() {
		IOUtil.assertFileIsReadable(INPUT_VCF);
		try (CloseableIterator<VariantContextDirectedBreakpoint> calls = getBreakpoints(INPUT_VCF)) {
			// Parallel populate cache requires
			// - indexed BCF/VCF OUTPUT.vcf.breakpoint.vcf
			// - (DONE) thread-safe cache 
			// Cache loading can then be performed in parallel per chunk
			populateCache_single_threaded(calls);
		}
	}
	private void populateCache_single_threaded(CloseableIterator<VariantContextDirectedBreakpoint> calls) {
		log.info("Loading variant evidence support");
		// TODO: resize this
		long readPairCount = getSamEvidenceSources().stream()
				.mapToLong(ses -> ses.getSVMetrics().STRUCTURAL_VARIANT_READ_PAIRS)
				.sum();
		long readCount = getSamEvidenceSources().stream()
				.mapToLong(ses -> ses.getSVMetrics().STRUCTURAL_VARIANT_READS)
				.sum();
		cache = new GreedyVariantAllocationCache(true, readPairCount, true, readCount, false, 0);
		try (CloseableIterator<DirectedEvidence> evidence = new AsyncBufferedIterator<>(getEvidenceIterator(), "mergedEvidence-cache")) {
			Iterator<VariantEvidenceSupport> annotator = new SequentialEvidenceAllocator(getContext(), calls, evidence, SAMEvidenceSource.maximumWindowSize(getContext(), getSamEvidenceSources(), getAssemblySource()), true);
			while (annotator.hasNext()) {
				VariantEvidenceSupport ves = annotator.next();
				for (DirectedEvidence e : ves.support) {
					if (e.isFromMultimappingFragment()) {
						cache.addBreakpoint((VariantContextDirectedBreakpoint)ves.variant, e);
					}
				}
			}
		}
	}
	private VariantContextDirectedBreakpoint annotate(VariantEvidenceSupport ves) {
		VariantCallingConfiguration vc = getContext().getConfig().getVariantCalling();
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), ves.variant);
		for (DirectedEvidence e : ves.support) {
			boolean shouldExclude = false;
			if (cache != null && e.isFromMultimappingFragment()) {
				shouldExclude = !cache.isBestBreakpoint((VariantContextDirectedBreakpoint)ves.variant, e);
			}
			if (!shouldExclude) {
				builder.addEvidence(e);
			}
		}
		VariantContextDirectedBreakpoint bp = (VariantContextDirectedBreakpoint)builder.make();
		if (!vc.writeFiltered) {
			if (bp.getBreakpointQual() < vc.minScore) return null;
			if (bp.getBreakpointEvidenceCount() < vc.minReads) return null;
			if (bp.isFiltered()) return null;
		}
		bp = (VariantContextDirectedBreakpoint)vc.applyConfidenceFilter(getContext(), bp);
		return bp;
	}
	public static void main(String[] argv) {
        System.exit(new AllocateEvidence().instanceMain(argv));
    }
}
