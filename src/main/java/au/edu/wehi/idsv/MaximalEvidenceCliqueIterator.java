package au.edu.wehi.idsv;

import java.util.Iterator;

import au.edu.wehi.idsv.graph.GraphNode;
import au.edu.wehi.idsv.graph.GraphNodeMergingIterator;
import au.edu.wehi.idsv.graph.RectangleGraphMaximalCliqueIterator;
import au.edu.wehi.idsv.graph.ScalingHelper;
import au.edu.wehi.idsv.util.WindowedSortingIterator;
import au.edu.wehi.idsv.vcf.VcfAttributes;
import au.edu.wehi.idsv.vcf.VcfSvConstants;

import com.google.common.base.Function;
import com.google.common.collect.AbstractIterator;
/**
 * Maximal clique summary evidence iterator 
 * @author Daniel Cameron
 *
 */
public class MaximalEvidenceCliqueIterator extends AbstractIterator<VariantContextDirectedEvidence> {
	public static final String MATE_BREAKEND_ID_SUFFIX_HIGH = "h";
	public static final String MATE_BREAKEND_ID_SUFFIX_LOW = "o";
	private VariantContextDirectedEvidence lastHigh = null;
	private final BreakendDirection targetLowDir;
	private final BreakendDirection targetHighDir;
	private final RectangleGraphMaximalCliqueIterator calc;
	private final ProcessingContext context;
	public MaximalEvidenceCliqueIterator(ProcessingContext processContext, Iterator<DirectedEvidence> evidenceIt, BreakendDirection lowDir, BreakendDirection highDir) {
		this.context = processContext;
		this.calc = new RectangleGraphMaximalCliqueIterator(
						// collapse evidence at the same location to a single node
						new GraphNodeMergingIterator(GraphNode.ByStartXY,
							// reorder due to soft clip margin changing record order
							new GraphNodeWindowedSortingIterator(context, 2 * processContext.getVariantCallingParameters().breakendMargin + 1,
								// convert evidence breakpoints to GraphNodes
								new EvidenceToGraphNodeIterator(evidenceIt))));
		this.targetLowDir = lowDir;
		this.targetHighDir = highDir;
	}
	private class GraphNodeWindowedSortingIterator extends WindowedSortingIterator<GraphNode> {
		public GraphNodeWindowedSortingIterator(final ProcessingContext processContext, final int windowSize, final Iterator<GraphNode> it) {
			super(it, new Function<GraphNode, Long>() {
				public Long apply(GraphNode arg) {
					return arg.startX;
				}
			}, windowSize, GraphNode.ByStartXY);
		}
	}
	private class EvidenceToGraphNodeIterator extends AbstractIterator<GraphNode> {
		private final Iterator<DirectedEvidence> it;
		public EvidenceToGraphNodeIterator(Iterator<DirectedEvidence> it) {
			this.it = it;
		}
		@Override
		protected GraphNode computeNext() {
			while (it.hasNext()) {
				DirectedEvidence evidence = it.next();
				GraphNode node = toGraphNode(evidence);
				if (node != null) {
					return node;
				}
			}
			return endOfData();
		}
	}
	private GraphNode toGraphNode(DirectedEvidence e) {
		BreakendSummary loc = e.getBreakendSummary();
		loc = context.getVariantCallingParameters().withMargin(context, loc);
		if (!(loc instanceof BreakpointSummary)) return null;
		
		BreakpointSummary bp = (BreakpointSummary)loc;
		assert(bp.referenceIndex >= 0);
		assert(bp.referenceIndex2 >= 0);
		assert(bp.start >= 1);
		assert(bp.start2 >= 1);
		long startX = context.getLinear().getLinearCoordinate(bp.referenceIndex, bp.start);
		long endX = startX + bp.end - bp.start;
		long startY = context.getLinear().getLinearCoordinate(bp.referenceIndex2, bp.start2);
		long endY = startY + bp.end2 - bp.start2;
		BreakendDirection lowDir = bp.direction;
		BreakendDirection highDir = bp.direction2;
		float weight = ((DirectedBreakpoint)e).getBreakpointQual();
		long scaledWeight = ScalingHelper.toScaledWeight(weight);
		if (scaledWeight <= 0) return null;
		GraphNode node = new GraphNode(startX, endX, startY, endY, scaledWeight);
		// Must have positive phred score  
		if (startX > startY) {
			// only take the lower half of the evidence since both sides of all breakpoints
			// have evidence
			// SC -> RemoteRealignedSoftClipEvidence
			// DP -> other half of the pair
			// Ass -> RemoteRealignedAssemblyEvidence
			return null;
			//node = node.flipAxis();
			//lowDir = bp.direction2;
			//highDir = bp.direction;
		}
		if (lowDir != targetLowDir || highDir != targetHighDir) return null;
		return node;
	}
	private VariantContextDirectedEvidence toVariant(GraphNode node, boolean isHighBreakend) {
		BreakpointSummary breakpoint = new BreakpointSummary(
				context.getLinear().getReferenceIndex(node.startX),
				targetLowDir,
				context.getLinear().getReferencePosition(node.startX),
				context.getLinear().getReferencePosition(node.endX),
				context.getLinear().getReferenceIndex(node.startY),
				targetHighDir,
				context.getLinear().getReferencePosition(node.startY),
				context.getLinear().getReferencePosition(node.endY));
		breakpoint = (BreakpointSummary)context.getVariantCallingParameters().withoutMargin(breakpoint);
		
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(context);
		// use a hash of the breakpoint as a (probably) unique identifier
		String id = String.format("gridss%d", Math.abs(breakpoint.hashCode()));
		builder.attribute(VcfSvConstants.BREAKEND_EVENT_ID_KEY, id);
		builder.attribute(VcfSvConstants.MATE_BREAKEND_ID_KEY, id + (isHighBreakend ? MATE_BREAKEND_ID_SUFFIX_LOW : MATE_BREAKEND_ID_SUFFIX_HIGH));
		builder.id(id + (isHighBreakend ? MATE_BREAKEND_ID_SUFFIX_HIGH : MATE_BREAKEND_ID_SUFFIX_LOW));
		if (isHighBreakend) {
			breakpoint = breakpoint.remoteBreakpoint();
		}
		builder.breakpoint(breakpoint, "");
		long scaledWeight = node.weight;
		double weight = ScalingHelper.toUnscaledWeight(scaledWeight);
		builder.phredScore(weight);
		builder.attribute(VcfAttributes.CALLED_QUAL, weight);
		VariantContextDirectedEvidence v = (VariantContextDirectedBreakpoint)builder.make();
		assert(v != null);
		return v;
	}
	@Override
	protected VariantContextDirectedEvidence computeNext() {
		if (lastHigh != null) {
			VariantContextDirectedEvidence result = lastHigh;
			lastHigh = null;
			return result;
		}
		if (calc.hasNext()) {
			GraphNode node = calc.next();
			VariantContextDirectedEvidence result = toVariant(node, false); 
			lastHigh = toVariant(node, true);
			return result;
		}
		return endOfData();
	}
}
