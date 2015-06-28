package au.edu.wehi.idsv.debruijn.positional;

import java.util.Iterator;
import java.util.NoSuchElementException;

import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.AssemblyParameters;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.SAMRecordAssemblyEvidence;

import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

/**
 * Assemblies non-reference breakend contigs
 * 
 * @author cameron.d
 *
 */
public class PositionalAssembler implements Iterator<SAMRecordAssemblyEvidence> {
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
	private void ensureAssembler() {
		while ((currentAssembler == null || !currentAssembler.hasNext()) && it.hasNext()) {
			// traverse contigs until we find one that has an assembly to call
			currentAssembler = createAssembler();
		}
	}
	private NonReferenceContigAssembler createAssembler() {
		AssemblyParameters ap = context.getAssemblyParameters();
		int maxEvidenceWidth = source.getMaxConcordantFragmentSize() - source.getMinConcordantFragmentSize() + 1;
		int maxReadLength = source.getMaxReadLength();
		int k = ap.k;
		int maxPathLength = ap.positionalMaxPathLengthInBases(maxEvidenceWidth);
		int maxPathCollapseLength = ap.positionalMaxPathCollapseLengthInBases(maxEvidenceWidth);
		int anchorAssemblyLength = ap.anchorAssemblyLength;
		int referenceIndex = it.peek().getBreakendSummary().referenceIndex;
		ReferenceIndexIterator evidenceIt = new ReferenceIndexIterator(it, referenceIndex);
		SupportNodeIterator supportIt = new SupportNodeIterator(k, evidenceIt, maxReadLength);
		EvidenceTracker trackedIt = new EvidenceTracker(supportIt);
		AggregateNodeIterator agIt = new AggregateNodeIterator(trackedIt);
		Iterator<KmerPathNode> pnIt = new PathNodeIterator(agIt, maxPathLength, k);
		if (ap.maxBaseMismatchForCollapse > 0) {
			pnIt = new PathCollapseIterator(pnIt, k, maxPathCollapseLength, ap.maxBaseMismatchForCollapse, ap.collapseBubblesOnly);
			pnIt = new PathSimplificationIterator(pnIt, maxPathLength, maxEvidenceWidth);
		}
		currentAssembler = new NonReferenceContigAssembler(pnIt, referenceIndex, maxEvidenceWidth, anchorAssemblyLength, k, source, trackedIt);
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
