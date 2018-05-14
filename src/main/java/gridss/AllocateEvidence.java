package gridss;

import java.util.Iterator;
import java.util.concurrent.ExecutorService;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;

import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.DirectedEvidenceOrder;
import au.edu.wehi.idsv.SAMEvidenceSource;
import au.edu.wehi.idsv.SequentialEvidenceAllocator;
import au.edu.wehi.idsv.SequentialEvidenceAllocator.VariantEvidenceSupport;
import au.edu.wehi.idsv.StructuralVariationCallBuilder;
import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;
import au.edu.wehi.idsv.VariantContextDirectedEvidence;
import au.edu.wehi.idsv.configuration.VariantCallingConfiguration;
import au.edu.wehi.idsv.util.AsyncBufferedIterator;
import au.edu.wehi.idsv.util.AutoClosingIterator;
import au.edu.wehi.idsv.validation.OrderAssertingIterator;
import au.edu.wehi.idsv.validation.PairedEvidenceTracker;
import gridss.cmdline.VcfTransformCommandLineProgram;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Log;

@CommandLineProgramProperties(
        summary = "Evidence reallocation is required to ensure that any given read/read pair/assembly, "
        		+ " supports only a consistent set breakpoints. \n"
        		+ "For discordant read pairs, this means that only a single breakpoint can be supported. \n"
        		+ "For split reads, only a single set of partial read mappings can be supported, and only one breakpoint per split can be supported.\n"
        		+ "For indels, only a single read mapping can be supported, and only one breakpoint per indel can be supported.\n"
        		+ "Note: the EID attribute must be populated with the relevant GRIDSS EvidenceID's using AllocateEvidence",  
        oneLineSummary = "Uniquely allocates evidence supporting multiple mutually exclusive breakpoints.",
        programGroup=gridss.cmdline.programgroups.VariantCalling.class
)
public class AllocateEvidence extends VcfTransformCommandLineProgram {
	private static final Log log = Log.getInstance(AllocateEvidence.class);
	@Argument(doc="Evidence allocation strategy used to uniquely assign evidence.")
	public EvidenceAllocationStrategy ALLOCATION_STRATEGY = EvidenceAllocationStrategy.GREEDY;
	public static enum EvidenceAllocationStrategy {
		GREEDY,
	}
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
	public CloseableIterator<VariantContextDirectedEvidence> iterator(CloseableIterator<VariantContextDirectedEvidence> calls, ExecutorService threadpool) {
		log.info("Allocating evidence"); 
		CloseableIterator<DirectedEvidence> evidence = new AsyncBufferedIterator<>(getEvidenceIterator(), "mergedEvidence-allocation");
		Iterator<VariantEvidenceSupport> annotator = new SequentialEvidenceAllocator(getContext(), calls, evidence, SAMEvidenceSource.maximumWindowSize(getContext(), getSamEvidenceSources(), getAssemblySource()), true);
		Iterator<VariantContextDirectedEvidence> it = Iterators.transform(annotator, bp -> annotate(bp));
		it = Iterators.filter(it, v -> v != null);
		return new AutoClosingIterator<>(it, calls, evidence);
	}
	private VariantContextDirectedEvidence annotate(VariantEvidenceSupport ves) {
		VariantCallingConfiguration vc = getContext().getConfig().getVariantCalling();
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), ves.variant);
		for (DirectedEvidence e : ves.support) {
			boolean shouldExclude = false;
			if (!shouldExclude) {
				builder.addEvidence(e);
			}
		}
		VariantContextDirectedEvidence be = builder.make();
		if (!vc.writeFiltered) {
			if (be.isFiltered()) return null;
			if (be.getPhredScaledQual() < vc.minScore) return null;
			if (be instanceof VariantContextDirectedBreakpoint) {
				VariantContextDirectedBreakpoint bp = (VariantContextDirectedBreakpoint)be;
				if (bp.getBreakpointReadCount() < vc.minReads) return null;
			} else {
				if (be.getBreakendReadCount() < vc.minReads) return null;
			}
		}
		be = vc.applyConfidenceFilter(getContext(), be);
		return be;
	}
	public static void main(String[] argv) {
        System.exit(new AllocateEvidence().instanceMain(argv));
    }
}
