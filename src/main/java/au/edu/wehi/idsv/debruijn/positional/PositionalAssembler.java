package au.edu.wehi.idsv.debruijn.positional;

import htsjdk.samtools.util.Log;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.NoSuchElementException;

import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.AssemblyParameters;
import au.edu.wehi.idsv.BreakendDirection;
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
	private final ProcessingContext context;
	private final AssemblyEvidenceSource source;
	private final PeekingIterator<DirectedEvidence> it;
	private final BreakendDirection direction;
	private NonReferenceContigAssembler currentAssembler = null;
	public PositionalAssembler(ProcessingContext context, AssemblyEvidenceSource source, Iterator<DirectedEvidence> backingIterator, BreakendDirection direction) {
		this.context = context;
		this.source = source;
		this.direction = direction;
		if (direction != null) {
			backingIterator = Iterators.filter(backingIterator, x -> x.getBreakendSummary() != null && x.getBreakendSummary().direction == this.direction);
		}
		this.it = Iterators.peekingIterator(backingIterator);
	}
	public PositionalAssembler(ProcessingContext context, AssemblyEvidenceSource source, Iterator<DirectedEvidence> it) {
		this(context, source, it, null);
	}
	public void exportGraph(File file) {
		if (currentAssembler != null) {
			currentAssembler.exportGraph(file);
		}
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
		int maxReadLength = source.getMaxMappedReadLength();
		int k = ap.k;
		int maxEvidenceDistance = maxSupportNodeWidth + maxReadLength + 2;
		int maxPathLength = ap.positionalMaxPathLengthInBases(maxReadLength);
		int maxPathCollapseLength = ap.positionalMaxPathCollapseLengthInBases(maxReadLength);
		int anchorAssemblyLength = ap.anchorAssemblyLength;
		int referenceIndex = it.peek().getBreakendSummary().referenceIndex;
		ReferenceIndexIterator evidenceIt = new ReferenceIndexIterator(it, referenceIndex);
		EvidenceTracker evidenceTracker = new EvidenceTracker();
		SupportNodeIterator supportIt = new SupportNodeIterator(k, evidenceIt, source.getMaxConcordantFragmentSize(), evidenceTracker);
		AggregateNodeIterator agIt = new AggregateNodeIterator(supportIt);
		Iterator<KmerNode> knIt = agIt;
		if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
			knIt = evidenceTracker.new AggregateNodeAssertionInterceptor(knIt);
		}
		PathNodeIterator pathNodeIt = new PathNodeIterator(knIt, maxPathLength, k); 
		Iterator<KmerPathNode> pnIt = pathNodeIt;
		if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
			pnIt = evidenceTracker.new PathNodeAssertionInterceptor(pnIt, "PathNodeIterator");
		}
		CollapseIterator collapseIt = null;
		PathSimplificationIterator simplifyIt = null;
		if (ap.maxBaseMismatchForCollapse > 0) {
			if (!ap.collapseBubblesOnly) {
				log.warn("Collapsing all paths is an exponential time operation. gridss is likely to hang if your genome contains repetative sequence");
				collapseIt = new PathCollapseIterator(pnIt, k, maxPathCollapseLength, ap.maxBaseMismatchForCollapse, false, 0);
			} else {
				collapseIt = new LeafBubbleCollapseIterator(pnIt, k, maxPathCollapseLength, ap.maxBaseMismatchForCollapse);
			}
			pnIt = collapseIt;
			if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
				pnIt = evidenceTracker.new PathNodeAssertionInterceptor(pnIt, "PathCollapseIterator");
			}
			simplifyIt = new PathSimplificationIterator(pnIt, maxPathLength, maxSupportNodeWidth);
			pnIt = simplifyIt;
			if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
				pnIt = evidenceTracker.new PathNodeAssertionInterceptor(pnIt, "PathSimplificationIterator");
			}
		}
		currentAssembler = new NonReferenceContigAssembler(pnIt, referenceIndex, maxEvidenceDistance, anchorAssemblyLength, k, source, evidenceTracker);
		if (ap.trackAlgorithmProgress && ap.debruijnGraphVisualisationDirectory != null) {
			ap.debruijnGraphVisualisationDirectory.mkdirs();
			String filename = String.format("positional-%s-%s.csv", context.getDictionary().getSequence(referenceIndex).getSequenceName(), direction);
			File file = new File(ap.debruijnGraphVisualisationDirectory, filename);
			PositionalDeBruijnGraphTracker exportTracker;
			try {
				exportTracker = new PositionalDeBruijnGraphTracker(file, supportIt, agIt, pathNodeIt, collapseIt, simplifyIt, evidenceTracker, currentAssembler);
				exportTracker.writeHeader();
				currentAssembler.setExportTracker(exportTracker);
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
