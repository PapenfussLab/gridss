package au.edu.wehi.idsv.debruijn.positional;

import htsjdk.samtools.util.Log;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.NoSuchElementException;

import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.AssemblyParameters;
import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.SAMRecordAssemblyEvidence;
import au.edu.wehi.idsv.visualisation.PositionalDeBruijnGraphTracker;

import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

/**
 * Assemblies non-reference breakend contigs
 * 
 * @author cameron.d
 *
 */
public class PositionalAssembler implements Iterator<SAMRecordAssemblyEvidence> {
	private static final Log log = Log.getInstance(PositionalAssembler.class);
	private ProcessingContext context;
	private AssemblyEvidenceSource source;
	private PeekingIterator<DirectedEvidence> it;
	private NonReferenceContigAssembler currentAssembler = null;
	public PositionalAssembler(ProcessingContext context, AssemblyEvidenceSource source, Iterator<DirectedEvidence> it) {
		this.context = context;
		this.source = source;
		this.it = Iterators.peekingIterator(it);
	}
	@Override
	public boolean hasNext() {
		ensureAssembler();
		return currentAssembler != null && currentAssembler.hasNext();
	}
	@Override
	public SAMRecordAssemblyEvidence next() {
		ensureAssembler();
		return currentAssembler.next();
	}
	private void flushIfRequired() {
		if (currentAssembler != null && !currentAssembler.hasNext()) {
			if (currentAssembler.getExportTracker() != null) {
				try {
					currentAssembler.getExportTracker().close();
				} catch (IOException e) {
					log.debug(e);
				}
			}
			currentAssembler = null;
		}
	}
	private void ensureAssembler() {
		flushIfRequired();
		while ((currentAssembler == null || !currentAssembler.hasNext()) && it.hasNext()) {
			// traverse contigs until we find one that has an assembly to call
			currentAssembler = createAssembler();
			flushIfRequired();
		}
	}
	private NonReferenceContigAssembler createAssembler() {
		AssemblyParameters ap = context.getAssemblyParameters();
		int maxSupportNodeWidth = source.getMaxConcordantFragmentSize() - source.getMinConcordantFragmentSize() + 1; 		
		int maxReadLength = source.getMaxReadLength();
		int k = ap.k;
		int maxEvidenceDistance = maxSupportNodeWidth + maxReadLength + 1;
		int maxPathLength = ap.positionalMaxPathLengthInBases(maxReadLength);
		int maxPathCollapseLength = ap.positionalMaxPathCollapseLengthInBases(maxReadLength);
		int anchorAssemblyLength = ap.anchorAssemblyLength;
		int referenceIndex = it.peek().getBreakendSummary().referenceIndex;
		ReferenceIndexIterator evidenceIt = new ReferenceIndexIterator(it, referenceIndex);
		SupportNodeIterator supportIt = new SupportNodeIterator(k, evidenceIt, source.getMaxConcordantFragmentSize());
		EvidenceTracker trackedIt = new EvidenceTracker(supportIt);
		AggregateNodeIterator agIt = new AggregateNodeIterator(trackedIt);
		Iterator<KmerNode> knIt = agIt;
		if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
			knIt = trackedIt.new AggregateNodeAssertionInterceptor(knIt);
		}
		PathNodeIterator pathNodeIt = new PathNodeIterator(knIt, maxPathLength, k); 
		Iterator<KmerPathNode> pnIt = pathNodeIt;
		if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
			pnIt = trackedIt.new PathNodeAssertionInterceptor(pnIt, "PathNodeIterator");
		}
		PathCollapseIterator collapseIt = null;
		PathSimplificationIterator simplifyIt = null;
		if (ap.maxBaseMismatchForCollapse > 0) {
			collapseIt = new PathCollapseIterator(pnIt, k, maxPathCollapseLength, ap.maxBaseMismatchForCollapse, ap.collapseBubblesOnly, 0); // TODO entropy parameter
			pnIt = collapseIt;
			if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
				pnIt = trackedIt.new PathNodeAssertionInterceptor(pnIt, "PathCollapseIterator");
			}
			simplifyIt = new PathSimplificationIterator(pnIt, maxPathLength, maxSupportNodeWidth);
			pnIt = simplifyIt;
			if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
				pnIt = trackedIt.new PathNodeAssertionInterceptor(pnIt, "PathSimplificationIterator");
			}
		}
		currentAssembler = new NonReferenceContigAssembler(pnIt, referenceIndex, maxEvidenceDistance, anchorAssemblyLength, k, source, trackedIt);
		if (ap.trackAlgorithmProgress && ap.debruijnGraphVisualisationDirectory != null) {
			ap.debruijnGraphVisualisationDirectory.mkdirs();
			String filename = String.format("positional-%s.csv", context.getDictionary().getSequence(referenceIndex).getSequenceName());
			File file = new File(ap.debruijnGraphVisualisationDirectory, filename);
			PositionalDeBruijnGraphTracker tracker;
			try {
				tracker = new PositionalDeBruijnGraphTracker(file, supportIt, agIt, pathNodeIt, collapseIt, simplifyIt, trackedIt, currentAssembler);
				tracker.writeHeader();
				currentAssembler.setExportTracker(tracker);
			} catch (IOException e) {
				log.debug(e);
			}
		}
		return currentAssembler;
	}
	private static class ReferenceIndexIterator implements PeekingIterator<DirectedEvidence> {
		private final PeekingIterator<DirectedEvidence> it;
		private final int referenceIndex;
		public ReferenceIndexIterator(PeekingIterator<DirectedEvidence> it, int referenceIndex) {
			this.it = it;
			this.referenceIndex = referenceIndex;
		}
		@Override
		public boolean hasNext() {
			return it.hasNext() && it.peek().getBreakendSummary().referenceIndex == referenceIndex;
		}

		@Override
		public DirectedEvidence next() {
			if (!hasNext()) throw new NoSuchElementException();
			return it.next();
		}

		@Override
		public DirectedEvidence peek() {
			if (hasNext()) return it.peek();
			throw new NoSuchElementException();
		}

		@Override
		public void remove() {
			throw new UnsupportedOperationException();
		}
		
	}
}
