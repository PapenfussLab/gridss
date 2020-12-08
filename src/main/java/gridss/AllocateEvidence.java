package gridss;

import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.*;
import au.edu.wehi.idsv.SequentialEvidenceAllocator.VariantEvidenceSupport;
import au.edu.wehi.idsv.configuration.VariantCallingConfiguration;
import au.edu.wehi.idsv.util.AsyncBufferedIterator;
import au.edu.wehi.idsv.util.AutoClosingIterator;
import au.edu.wehi.idsv.util.AutoClosingMergedIterator;
import au.edu.wehi.idsv.validation.OrderAssertingIterator;
import au.edu.wehi.idsv.validation.PairedEvidenceTracker;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;
import gridss.cmdline.VcfTransformCommandLineProgram;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;

import java.io.Closeable;
import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.ExecutorService;

@CommandLineProgramProperties(
        summary = "Evidence reallocation is required to ensure that any given read/read pair/assembly "
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
	@Argument(doc="Indicates whether supporting assemblies should be allocated.")
	public boolean ALLOCATE_ASSEMBLIES = true;
	@Argument(doc="Indicates whether supporting reads should be allocated.")
	public boolean ALLOCATE_READS = true;
	public static enum EvidenceAllocationStrategy {
		GREEDY,
	}
	private CalledBreakpointPositionLookup lookup = new CalledBreakpointPositionLookup();
	public CloseableIterator<DirectedEvidence> getReadIterator() {
		CloseableIterator<DirectedEvidence> evidenceIt;
		List<SAMEvidenceSource> sources = getSamEvidenceSources();
		sources.stream().forEach(ses -> ses.assertPreprocessingComplete());
		evidenceIt = SAMEvidenceSource.mergedIterator(ImmutableList.<SAMEvidenceSource>builder().addAll(sources).build(), true, SAMEvidenceSource.EvidenceSortOrder.EvidenceStartPosition);
		if (Defaults.SANITY_CHECK_ITERATORS) {
			evidenceIt = new AutoClosingIterator<>(
					new PairedEvidenceTracker<>("Reads", new OrderAssertingIterator<>(evidenceIt, DirectedEvidenceOrder.ByNatural), false),
					evidenceIt);
		}
		return evidenceIt;
	}
	public CloseableIterator<DirectedEvidence> getAssemblyIterator() {
		CloseableIterator<DirectedEvidence> evidenceIt;
		evidenceIt = new AggregateEvidenceSource(getContext(), getAssemblySource(), null, SAMEvidenceSource.EvidenceSortOrder.EvidenceStartPosition).iterator();
		if (Defaults.SANITY_CHECK_ITERATORS) {
			evidenceIt = new AutoClosingIterator<>(
					new PairedEvidenceTracker<>("Assemblies", new OrderAssertingIterator<>(evidenceIt, DirectedEvidenceOrder.ByNatural), false),
					evidenceIt);
		}
		return evidenceIt;
	}
	@Override
	public CloseableIterator<VariantContextDirectedEvidence> iterator(CloseableIterator<VariantContextDirectedEvidence> calls, ExecutorService threadpool) {
		log.info("Allocating evidence"); 
		CloseableIterator<DirectedEvidence> rawReads = new AsyncBufferedIterator<>(getReadIterator(), "mergedReads-allocation");
		CloseableIterator<DirectedEvidence> reads = new AsyncBufferedIterator<>(annotateAssembly(rawReads), "annotate-associated-assembly");
		CloseableIterator<DirectedEvidence> assemblies = new AsyncBufferedIterator<>(getAssemblyIterator(), "assembly-allocation");
		Iterator<VariantEvidenceSupport> annotator = new SequentialEvidenceAllocator(getContext(), calls, reads, assemblies, SAMEvidenceSource.maximumWindowSize(getContext(), getSamEvidenceSources(), getAssemblySource()), true);
		CloseableIterator<VariantEvidenceSupport> bufferedAnnotator = new AsyncBufferedIterator<>(annotator, "annotator", 2, 8);
		Iterator<VariantContextDirectedEvidence> it = Iterators.transform(bufferedAnnotator, bp -> annotate(bp));
		it = Iterators.filter(it, v -> v != null);
		return new AutoClosingIterator<>(it, calls, rawReads, reads, assemblies, bufferedAnnotator);
	}
	private CloseableIterator<DirectedEvidence> annotateAssembly(CloseableIterator<DirectedEvidence> it) {
		List<Closeable> assToClose = new ArrayList<>();
		List<Iterator<SAMRecord>> rawAssemblies = new ArrayList<>();
		int windowSize = 0;
		for (AssemblyEvidenceSource aes : getAssemblySource()) {
			// need to use the raw breakend assembly file (prior to realignment) so we annotate correctly
			File assemblyFile = aes.getFile();
			if (assemblyFile == null || !assemblyFile.exists()) {
				if (getContext().getConfig().getVariantCalling().ignoreMissingAssemblyFile) {
					log.error("Missing assembly file. BAN* annotations will be incorrect.");
					return it;
				} else {
					throw new RuntimeException("Missing assembly file " + (assemblyFile == null ? "" : assemblyFile.getName()));
				}
			}
			windowSize = Math.max(windowSize, aes.getMaxAssemblyLength() + 2 * aes.getMaxConcordantFragmentSize());
			// defensive over-eager loading
			windowSize *= 2;
			SamReader reader = getContext().getSamReader(assemblyFile);
			SAMRecordIterator assit = reader.iterator();
			rawAssemblies.add(assit);
			assToClose.add(assit);
			assToClose.add(reader);
		}
		AutoClosingMergedIterator mergedAssemblies = new AutoClosingMergedIterator(rawAssemblies, new SAMRecordCoordinateOnlyComparator());
		return new AutoClosingIterator<>(new AssemblyAssociator(it, mergedAssemblies, windowSize), assToClose.toArray(new Closeable[0]));
	}
	private VariantContextDirectedEvidence annotate(VariantEvidenceSupport ves) {
		VariantCallingConfiguration vc = getContext().getConfig().getVariantCalling();
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), lookup, ves.variant);
		builder.setUpdateAssemblyInformation(ALLOCATE_ASSEMBLIES);
		builder.setUpdateReadInformation(ALLOCATE_READS);
		for (DirectedEvidence e : ves.support) {
			builder.addEvidence(e);
		}
		VariantContextDirectedEvidence be = builder.make();
		if (ALLOCATE_READS && ALLOCATE_ASSEMBLIES) {
			if (!vc.writeFiltered) {
				if (be.isFiltered()) return null;
				if (be.getPhredScaledQual() < vc.minScore) return null;
				if (be instanceof VariantContextDirectedBreakpoint) {
					VariantContextDirectedBreakpoint bp = (VariantContextDirectedBreakpoint) be;
					if (bp.getBreakpointSupportingFragmentCount() < vc.minReads) return null;
				} else {
					if (be.getBreakendSupportingFragmentCount() < vc.minReads) return null;
				}
			}
			be = vc.applyConfidenceFilter(getContext(), be);
		}
		return be;
	}
	public static void main(String[] argv) {
        System.exit(new AllocateEvidence().instanceMain(argv));
    }
}
