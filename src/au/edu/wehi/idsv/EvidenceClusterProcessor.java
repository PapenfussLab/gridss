package au.edu.wehi.idsv;

import java.util.Iterator;
import java.util.PriorityQueue;

import au.edu.wehi.idsv.graph.GraphNode;
import au.edu.wehi.idsv.graph.MaximalClique;
import au.edu.wehi.idsv.vcf.VcfAttributes;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;
import com.google.common.collect.UnmodifiableIterator;
/**
 * Calls breakpoints from the given evidence
 * 
 * @author Daniel Cameron
 */
public class EvidenceClusterProcessor implements Iterable<VariantContextDirectedEvidence> {
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
			GraphNode node = new GraphNode(startX, endX, startY, endY, (float)Models.llr(evidence));
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
	public Iterator<VariantContextDirectedEvidence> iterator() {
		UnmodifiableIterator<VariantContextDirectedEvidence> x = Iterators.mergeSorted(ImmutableList.of(
				new EvidenceClusterProcessorMaximalCliqueIterator(ff, BreakendDirection.Forward, BreakendDirection.Forward),
				new EvidenceClusterProcessorMaximalCliqueIterator(fb, BreakendDirection.Forward, BreakendDirection.Backward),
				new EvidenceClusterProcessorMaximalCliqueIterator(bf, BreakendDirection.Backward, BreakendDirection.Forward),
				new EvidenceClusterProcessorMaximalCliqueIterator(bb, BreakendDirection.Backward, BreakendDirection.Backward)),
			IdsvVariantContext.ByLocationStart);
		// TODO: filtering and processing of the maximal clique
		// don't call weak evidence
		// don't call if there is a much stronger call nearby
		return x;
	}
	/**
	 * Maximal clique summary evidence iterator 
	 * @author Daniel Cameron
	 *
	 */
	public class EvidenceClusterProcessorMaximalCliqueIterator extends AbstractIterator<VariantContextDirectedEvidence> { 
		private final PriorityQueue<GraphNode> highBreakend = new PriorityQueue<GraphNode>(GraphNode.ByStartYXEndYX);
		private final PeekingIterator<GraphNode> it;
		private final BreakendDirection lowDir;
		private final BreakendDirection highDir;
		public EvidenceClusterProcessorMaximalCliqueIterator(MaximalClique completeGraph, BreakendDirection lowDir, BreakendDirection highDir) {
			this.it = Iterators.peekingIterator(completeGraph.getAllMaximalCliques());
			this.lowDir = lowDir;
			this.highDir = highDir;
		}
		private VariantContextDirectedEvidence toVariant(GraphNode node, boolean isHighBreakend) {
			BreakpointSummary breakpoint = new BreakpointSummary(
					context.getLinear().getReferenceIndex(node.startX),
					lowDir,
					context.getLinear().getReferencePosition(node.startX),
					context.getLinear().getReferencePosition(node.endX),
					context.getLinear().getReferenceIndex(node.startY),
					highDir,
					context.getLinear().getReferencePosition(node.startY),
					context.getLinear().getReferencePosition(node.endY));
			IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(context);
			// use a hash of the breakpoint as a (probably) unique identifier 
			builder.id(String.format("idsv%d", Math.abs(breakpoint.hashCode())));
			if (isHighBreakend) {
				breakpoint = breakpoint.remoteBreakpoint();
			}
			builder.breakpoint(breakpoint, null);
			builder.phredScore(node.weight);
			builder.attribute(VcfAttributes.LOG_LIKELIHOOD_RATIO_BREAKPOINT, node.weight);
			return (VariantContextDirectedEvidence)builder.make();
		}
		@Override
		protected VariantContextDirectedEvidence computeNext() {
			if (!it.hasNext()) {
				if (highBreakend.isEmpty()) {
					return endOfData();
				} else {
					return toVariant(highBreakend.poll(), true);
				}
			}
			// BUG: sort order of output is by location end, not start
			// calls for breakend mapping to two difference location can be returned
			// out of order. To fix, we need to know the max call window size (=max frag size)
			// and buffer the input iterator by that much to swap the sort order
			// Alteratively, MaximalClique can traverse backwards which would return cliques
			// in starting order, instead of ending order. Easiest way to do this is to make
			// all coordinates given to MaximalClique negative
			if (highBreakend.isEmpty() || highBreakend.peek().startY > it.peek().startX) {
				// iterator is next
				GraphNode node = it.next();
				highBreakend.add(node);
				return toVariant(node, false);
			} else {
				return toVariant(highBreakend.poll(), true);
			}
		}
	}
}
 