package au.edu.wehi.socrates;

import java.util.Iterator;

import au.edu.wehi.socrates.graph.TrapezoidGraph;
import au.edu.wehi.socrates.graph.TrapezoidGraphNode;

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
public class EvidenceClusterProcessor implements Iterable<BreakpointLocation> {
	private final TrapezoidGraph ff = new TrapezoidGraph();
	private final TrapezoidGraph fb = new TrapezoidGraph();
	private final TrapezoidGraph bf = new TrapezoidGraph();
	private final TrapezoidGraph bb = new TrapezoidGraph();
	private final ProcessingContext context;
	public EvidenceClusterProcessor(ProcessingContext context) {
		this.context = context;
	}
	public void addEvidence(DirectedEvidence evidence) {
		BreakpointLocation loc = evidence.getBreakpointLocation();
		if (evidence instanceof SoftClipEvidence) {
			// TODO: expand evidence region for soft-clips
			// Socrates expands by +- 3bp
		}
		if (loc instanceof BreakpointInterval) {
			BreakpointInterval interval = (BreakpointInterval)loc;
			long startX = context.getLinear().getLinearCoordinate(interval.referenceIndex, interval.start);
			long endX = startX + interval.end - interval.start;
			long startY = context.getLinear().getLinearCoordinate(interval.referenceIndex2, interval.start2);
			long endY = startY + interval.end2 - interval.start2;
			BreakpointDirection lowDir = interval.direction;
			BreakpointDirection highDir = interval.direction2;
			if (startX > startY) {
				// reverse evidence: need to flip coordinates
				long tmp = startX;
				startX = startY;
				startY = tmp;
				tmp = endX;
				endX = endY;
				endY = tmp;
				BreakpointDirection tmpDir = lowDir;
				lowDir = highDir;
				highDir = tmpDir;
			}
			TrapezoidGraphNode node = new TrapezoidGraphNode(startX, endX, startY, endY, loc.qual);
			if (lowDir == BreakpointDirection.Forward) {
				if (highDir == BreakpointDirection.Forward) {
					ff.add(node);
				} else {
					fb.add(node);
				}
			} else {
				if (highDir == BreakpointDirection.Forward) {
					bf.add(node);
				} else {
					bb.add(node);
				}
			}
		} else {
			// TODO: process evidence such as short SCs and OEAs which do not map to a remote location
		}
	}
	@Override
	public Iterator<BreakpointLocation> iterator() {
		UnmodifiableIterator<BreakpointInterval> x = Iterators.mergeSorted(ImmutableList.of(
				new EvidenceClusterProcessorMaximalCliqueIterator(ff, BreakpointDirection.Forward, BreakpointDirection.Forward),
				new EvidenceClusterProcessorMaximalCliqueIterator(fb, BreakpointDirection.Forward, BreakpointDirection.Backward),
				new EvidenceClusterProcessorMaximalCliqueIterator(bf, BreakpointDirection.Backward, BreakpointDirection.Forward),
				new EvidenceClusterProcessorMaximalCliqueIterator(bb, BreakpointDirection.Backward, BreakpointDirection.Backward)),
				new Ordering<BreakpointInterval>() {
					public int compare(BreakpointInterval o1, BreakpointInterval o2) {
						// same order as our source iterators
						return ComparisonChain.start()
								.compare(o1.end, o2.end)
								.compare(o1.start2, o2.start2)
								.result();
					}
		});
		// downcast iterator to BreakpointLocation
		return Iterators.<BreakpointLocation>concat(x);
		// TODO: filtering and processing of the maximal clique
		// don't call weak evidence
		// don't call if there is a much stronger call nearby
	}
	public class EvidenceClusterProcessorMaximalCliqueIterator implements Iterator<BreakpointInterval> {
		private final Iterator<TrapezoidGraphNode> it;
		private final BreakpointDirection lowDir;
		private final BreakpointDirection highDir;
		public EvidenceClusterProcessorMaximalCliqueIterator(TrapezoidGraph completeGraph, BreakpointDirection lowDir, BreakpointDirection highDir) {
			this.it = completeGraph.getAllMaximalCliques();
			this.lowDir = lowDir;
			this.highDir = highDir;
		}
		@Override
		public boolean hasNext() {
			return it.hasNext();
		}
		@Override
		public BreakpointInterval next() {
			TrapezoidGraphNode node = it.next();
			return new BreakpointInterval(
					context.getLinear().getReferenceIndex(node.startX),
					lowDir,
					context.getLinear().getReferencePosition(node.startX),
					context.getLinear().getReferencePosition(node.endX),
					context.getLinear().getReferenceIndex(node.startY),
					highDir,
					context.getLinear().getReferencePosition(node.startY),
					context.getLinear().getReferencePosition(node.endY),
					node.weight);
		}
		@Override
		public void remove() {
			it.remove();
		}
	}
}
 