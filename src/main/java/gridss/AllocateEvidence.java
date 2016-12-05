package gridss;

import java.util.Iterator;
import java.util.concurrent.ExecutorService;

import com.google.common.collect.Iterators;

import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.GreedyVariantAllocationCache;
import au.edu.wehi.idsv.SAMEvidenceSource;
import au.edu.wehi.idsv.SequentialEvidenceAllocator;
import au.edu.wehi.idsv.SequentialEvidenceAllocator.VariantEvidenceSupport;
import au.edu.wehi.idsv.StructuralVariationCallBuilder;
import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;
import au.edu.wehi.idsv.configuration.VariantCallingConfiguration;
import au.edu.wehi.idsv.util.AsyncBufferedIterator;
import au.edu.wehi.idsv.util.AutoClosingIterator;
import gridss.cmdline.VcfTransformCommandLineProgram;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;

@CommandLineProgramProperties(
        usage = "Evidence reallocation is required to ensure that any given read/read pair/assembly, "
        		+ " supports only a consistent set breakpoints. \n"
        		+ "For discordant read pairs, this means that only a single breakpoint can be supported. \n"
        		+ "For split reads, only a single set of partial read mappings can be supported, and only one breakpoint per split can be supported.\n"
        		+ "For indels, only a single read mapping can be supported, and only one breakpoint per indel can be supported.\n"
        		+ "Note: the EID attribute must be populated with the relevant GRIDSS EvidenceID's using AllocateEvidence",  
        usageShort = "Unique allocates evidence supporting multiple mutually exclusive breakpoints."
)
public class AllocateEvidence extends VcfTransformCommandLineProgram {
	private static final Log log = Log.getInstance(AllocateEvidence.class);
	/**
	 * Greedily assigning each evidence to the best variant is fine since
	 * we don't need to consider sub-optimal locations.
	 * This reduces the size of the lookup cache required
	 */
	private static final boolean ALLOCATE_TO_BEST = true;
	
	@Option(doc="Evidence allocation strategy used to uniquely assign evidence.")
	public EvidenceAllocationStrategy ALLOCATION_STRATEGY = EvidenceAllocationStrategy.GREEDY;
	public static enum EvidenceAllocationStrategy {
		GREEDY,
	}
	private GreedyVariantAllocationCache cache;
	@Override
	public CloseableIterator<VariantContextDirectedBreakpoint> iterator(CloseableIterator<VariantContextDirectedBreakpoint> calls, ExecutorService threadpool) {
		populateCache();
		log.info("Allocating evidence");
		CloseableIterator<DirectedEvidence> evidence = new AsyncBufferedIterator<>(getEvidenceIterator(), "mergedEvidence-allocation");
		Iterator<VariantEvidenceSupport> annotator = new SequentialEvidenceAllocator(getContext(), calls, evidence, SAMEvidenceSource.maximumWindowSize(getContext(), getSamEvidenceSources(), getAssemblySource()), ALLOCATE_TO_BEST);
		Iterator<VariantContextDirectedBreakpoint> it = Iterators.transform(annotator, bp -> annotate(bp));
		it = Iterators.filter(it, v -> v != null);
		return new AutoClosingIterator<>(it, calls, evidence);
	}
	public void populateCache() {
		IOUtil.assertFileIsReadable(INPUT_VCF);
		try (CloseableIterator<VariantContextDirectedBreakpoint> calls = getBreakpoints(INPUT_VCF)) {
			populateCache(calls);
		}
	}
	public void populateCache(CloseableIterator<VariantContextDirectedBreakpoint> calls) {
		log.info("Loading variant evidence support");
		cache = new GreedyVariantAllocationCache();
		try (CloseableIterator<DirectedEvidence> evidence = new AsyncBufferedIterator<>(getEvidenceIterator(), "mergedEvidence-cache")) {
			Iterator<VariantEvidenceSupport> annotator = new SequentialEvidenceAllocator(getContext(), calls, evidence, SAMEvidenceSource.maximumWindowSize(getContext(), getSamEvidenceSources(), getAssemblySource()), ALLOCATE_TO_BEST);
			while (annotator.hasNext()) {
				VariantEvidenceSupport ves = annotator.next();
				for (DirectedEvidence e : ves.support) {
					cache.addBreakpoint((VariantContextDirectedBreakpoint)ves.variant, e);
				}
			}
		}
	}
	private VariantContextDirectedBreakpoint annotate(VariantEvidenceSupport ves) {
		VariantCallingConfiguration vc = getContext().getConfig().getVariantCalling();
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), ves.variant);
		for (DirectedEvidence e : ves.support) {
			if (cache.isBestBreakpoint((VariantContextDirectedBreakpoint)ves.variant, e)) {
				builder.addEvidence(e);
			}
		}
		VariantContextDirectedBreakpoint bp = (VariantContextDirectedBreakpoint)builder.make();
		if (!vc.writeFiltered) {
			if (bp.getBreakpointQual() < vc.minScore) return null;
			if (bp.isFiltered()) return null;
		}
		bp = (VariantContextDirectedBreakpoint)vc.applyConfidenceFilter(getContext(), bp);
		return bp;
	}
}
