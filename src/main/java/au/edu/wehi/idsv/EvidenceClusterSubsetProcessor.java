package au.edu.wehi.idsv;

import java.util.Iterator;

import com.google.common.collect.AbstractIterator;
/**
 * Calls breakpoints from the given evidence
 * 
 * @author Daniel Cameron
 */
public class EvidenceClusterSubsetProcessor extends EvidenceClusterProcessor {
	public EvidenceClusterSubsetProcessor(
			ProcessingContext context,
			Iterator<DirectedEvidence> evidenceIt,
			int fromReferenceIndex,
			int toReferenceIndex) {
		super(context, new ChromosomeFilteringIterator(evidenceIt, Math.min(fromReferenceIndex, toReferenceIndex), Math.max(fromReferenceIndex, toReferenceIndex)));
	}
	protected static class ChromosomeFilteringIterator extends AbstractIterator<DirectedEvidence> {
		private int fromReferenceIndex;
		private int toReferenceIndex;
		private Iterator<DirectedEvidence> it;
		public ChromosomeFilteringIterator(
				Iterator<DirectedEvidence> evidenceIt,
				int fromReferenceIndex,
				int toReferenceIndex) {
			this.it = evidenceIt;
			this.fromReferenceIndex = fromReferenceIndex;
			this.toReferenceIndex = toReferenceIndex;
		}
		@Override
		protected DirectedEvidence computeNext() {
			while (it.hasNext()) {
				DirectedEvidence e = it.next();
				if (!isFiltered(e)) return e;
			}
			return endOfData();
		}
		protected boolean isFiltered(DirectedEvidence e) {
			BreakendSummary loc = e.getBreakendSummary();
			if (loc instanceof BreakpointSummary) {
				return isFiltered((BreakpointSummary)loc);
			}
			return isFiltered(loc);
		}
		protected boolean isFiltered(BreakendSummary loc) {
			return loc.referenceIndex != fromReferenceIndex && loc.referenceIndex != toReferenceIndex;
		}
		protected boolean isFiltered(BreakpointSummary loc) {
			return !(loc.referenceIndex == fromReferenceIndex && loc.referenceIndex2 == toReferenceIndex)
					&& !(loc.referenceIndex == toReferenceIndex && loc.referenceIndex2 == fromReferenceIndex);
		}
	}
}
 