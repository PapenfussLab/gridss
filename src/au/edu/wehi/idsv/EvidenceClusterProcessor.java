package au.edu.wehi.idsv;

import java.util.Iterator;

import au.edu.wehi.idsv.graph.GraphNode;
import au.edu.wehi.idsv.graph.MaximalClique;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;
import com.google.common.collect.Ordering;
import com.google.common.collect.UnmodifiableIterator;
/**
 * Calls breakpoints from the given evidence
 * 
 * @author Daniel Cameron
 */
public class EvidenceClusterProcessor implements Iterable<BreakpointSummary> {
	private final MaximalClique ff = new MaximalClique();
	private final MaximalClique fb = new MaximalClique();
	private final MaximalClique bf = new MaximalClique();
	private final MaximalClique bb = new MaximalClique();
	private final ProcessingContext context;
	public EvidenceClusterProcessor(ProcessingContext context) {
		this.context = context;
	}
	public void addEvidence(DirectedEvidence evidence) {
		BreakendSummary loc = evidence.getBreakendSummary();
		if (evidence instanceof SoftClipEvidence) {
			// TODO: expand evidence region for soft-clips
			// Socrates expands by +- 3bp
		}
		if (loc instanceof BreakpointSummary) {
			BreakpointSummary interval = (BreakpointSummary)loc;
			if (filterOut(interval)) return;
			long startX = context.getLinear().getLinearCoordinate(interval.referenceIndex, interval.start);
			long endX = startX + interval.end - interval.start;
			long startY = context.getLinear().getLinearCoordinate(interval.referenceIndex2, interval.start2);
			long endY = startY + interval.end2 - interval.start2;
			BreakendDirection lowDir = interval.direction;
			BreakendDirection highDir = interval.direction2;
			if (startX > startY) {
				// reverse evidence: need to flip coordinates
				long tmp = startX;
				startX = startY;
				startY = tmp;
				tmp = endX;
				endX = endY;
				endY = tmp;
				BreakendDirection tmpDir = lowDir;
				lowDir = highDir;
				highDir = tmpDir;
			}
			GraphNode node = new GraphNode(startX, endX, startY, endY, (float)loc.evidence.getScore());
			addNode(lowDir, highDir, node);
		} else {
			if (filterOut(loc)) return;
			// TODO: process evidence such as short SCs and OEAs which do not map to a remote location
		}
	}
	/**
	 * Determine whether the given evidence should be filtered
	 * @param loc evidence location to consider filtering
	 * @return true if the evidence should be filtered out 
	 */
	protected boolean filterOut(BreakendSummary loc) {
		return false;
	}
	/**
	 * Determine whether the given evidence should be filtered
	 * @param loc evidence location to consider filtering
	 * @return true if the evidence should be filtered out 
	 */
	protected boolean filterOut(BreakpointSummary loc) {
		return false;
	}
	private void addNode(BreakendDirection lowDir, BreakendDirection highDir, GraphNode node) {
		if (lowDir == BreakendDirection.Forward) {
			if (highDir == BreakendDirection.Forward) {
				ff.add(node);
			} else {
				fb.add(node);
			}
		} else {
			if (highDir == BreakendDirection.Forward) {
				bf.add(node);
			} else {
				bb.add(node);
			}
		}
	}
	@Override
	public Iterator<BreakpointSummary> iterator() {
		UnmodifiableIterator<BreakpointSummary> x = Iterators.mergeSorted(ImmutableList.of(
				new EvidenceClusterProcessorMaximalCliqueIterator(ff, BreakendDirection.Forward, BreakendDirection.Forward),
				new EvidenceClusterProcessorMaximalCliqueIterator(fb, BreakendDirection.Forward, BreakendDirection.Backward),
				new EvidenceClusterProcessorMaximalCliqueIterator(bf, BreakendDirection.Backward, BreakendDirection.Forward),
				new EvidenceClusterProcessorMaximalCliqueIterator(bb, BreakendDirection.Backward, BreakendDirection.Backward)),
				new Ordering<BreakpointSummary>() {
					public int compare(BreakpointSummary o1, BreakpointSummary o2) {
						// same order as our source iterators
						return ComparisonChain.start()
								.compare(o1.end, o2.end)
								.compare(o1.start2, o2.start2)
								.result();
					}
		});
		return x;
		// TODO: filtering and processing of the maximal clique
		// don't call weak evidence
		// don't call if there is a much stronger call nearby
	}
	public class EvidenceClusterProcessorMaximalCliqueIterator implements Iterator<BreakpointSummary> {
		private final Iterator<GraphNode> it;
		private final BreakendDirection lowDir;
		private final BreakendDirection highDir;
		public EvidenceClusterProcessorMaximalCliqueIterator(MaximalClique completeGraph, BreakendDirection lowDir, BreakendDirection highDir) {
			this.it = completeGraph.getAllMaximalCliques();
			this.lowDir = lowDir;
			this.highDir = highDir;
		}
		@Override
		public boolean hasNext() {
			return it.hasNext();
		}
		@Override
		public BreakpointSummary next() {
			GraphNode node = it.next();
			return new BreakpointSummary(
					context.getLinear().getReferenceIndex(node.startX),
					lowDir,
					context.getLinear().getReferencePosition(node.startX),
					context.getLinear().getReferencePosition(node.endX),
					context.getLinear().getReferenceIndex(node.startY),
					highDir,
					context.getLinear().getReferencePosition(node.startY),
					context.getLinear().getReferencePosition(node.endY),
					new EvidenceMetrics(node.evidence));
		}
		@Override
		public void remove() {
			it.remove();
		}
	}
}
 