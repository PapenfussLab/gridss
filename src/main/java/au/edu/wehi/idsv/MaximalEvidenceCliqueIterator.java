package au.edu.wehi.idsv;

import java.util.Iterator;

import au.edu.wehi.idsv.graph.RectangleGraphMaximalCliqueIterator;
import au.edu.wehi.idsv.graph.RectangleGraphNode;
import au.edu.wehi.idsv.graph.RectangleGraphNodeMergingIterator;
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
	private final VariantIdGenerator idGenerator;
	public MaximalEvidenceCliqueIterator(ProcessingContext processContext, Iterator<DirectedEvidence> evidenceIt, BreakendDirection lowDir, BreakendDirection highDir, VariantIdGenerator idGenerator) {
		this.context = processContext;
		this.calc = new RectangleGraphMaximalCliqueIterator(
						// collapse evidence at the same location to a single node
						new RectangleGraphNodeMergingIterator(RectangleGraphNode.ByStartXY,
							// reorder due to soft clip margin changing record order
							new GraphNodeWindowedSortingIterator(context, 2 * processContext.getVariantCallingParameters().breakendMargin + 1,
								// convert evidence breakpoints to GraphNodes
								new EvidenceToGraphNodeIterator(evidenceIt))));
		this.targetLowDir = lowDir;
		this.targetHighDir = highDir;
		this.idGenerator = idGenerator;
	}
	private class GraphNodeWindowedSortingIterator extends WindowedSortingIterator<RectangleGraphNode> {
		public GraphNodeWindowedSortingIterator(final ProcessingContext processContext, final int windowSize, final Iterator<RectangleGraphNode> it) {
			super(it, new Function<RectangleGraphNode, Long>() {
				public Long apply(RectangleGraphNode arg) {
					return arg.startX;
				}
			}, windowSize, RectangleGraphNode.ByStartXY);
		}
	}
	private class EvidenceToGraphNodeIterator extends AbstractIterator<RectangleGraphNode> {
		private final Iterator<DirectedEvidence> it;
		public EvidenceToGraphNodeIterator(Iterator<DirectedEvidence> it) {
			this.it = it;
		}
		@Override
		protected RectangleGraphNode computeNext() {
			while (it.hasNext()) {
				DirectedEvidence evidence = it.next();
				RectangleGraphNode node = toGraphNode(evidence);
				if (node != null) {
					return node;
				}
			}
			return endOfData();
		}
	}
	private RectangleGraphNode toGraphNode(DirectedEvidence e) {
		BreakendSummary loc = e.getBreakendSummary();
		loc = context.getVariantCallingParameters().withMargin(context, loc);
		if (!(loc instanceof BreakpointSummary)) return null;		
		BreakpointSummary bp = (BreakpointSummary)loc;
		long startX = context.getLinear().getLinearCoordinate(bp.referenceIndex, bp.start);
		long endX = startX + bp.end - bp.start;
		long startY = context.getLinear().getLinearCoordinate(bp.referenceIndex2, bp.start2);
		long endY = startY + bp.end2 - bp.start2;
		BreakendDirection lowDir = bp.direction;
		BreakendDirection highDir = bp.direction2;
		float weight = ((DirectedBreakpoint)e).getBreakpointQual();
		long scaledWeight = ScalingHelper.toScaledWeight(weight);
		if (scaledWeight <= 0) return null;
		RectangleGraphNode node = new RectangleGraphNode(startX, endX, startY, endY, scaledWeight);
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
	private VariantContextDirectedEvidence toVariant(String event, RectangleGraphNode node, BreakpointSummary breakpoint, boolean isHighBreakend) {
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(context);
		builder.attribute(VcfSvConstants.BREAKEND_EVENT_ID_KEY, event);
		builder.attribute(VcfSvConstants.MATE_BREAKEND_ID_KEY, event + (isHighBreakend ? MATE_BREAKEND_ID_SUFFIX_LOW : MATE_BREAKEND_ID_SUFFIX_HIGH));
		builder.id(event + (isHighBreakend ? MATE_BREAKEND_ID_SUFFIX_HIGH : MATE_BREAKEND_ID_SUFFIX_LOW));
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
	private BreakpointSummary toBreakpointSummary(RectangleGraphNode node) {
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
		// sanity check that the resultant breakpoint makes sense
		assert(breakpoint.start >= 1);
		assert(breakpoint.start2 >= 1);
		assert(breakpoint.referenceIndex >= 0);
		assert(breakpoint.referenceIndex2 >= 0);
		return breakpoint;
	}
	@Override
	protected VariantContextDirectedEvidence computeNext() {
		if (lastHigh != null) {
			VariantContextDirectedEvidence result = lastHigh;
			lastHigh = null;
			return result;
		}
		if (calc.hasNext()) {
			RectangleGraphNode node = calc.next();
			BreakpointSummary breakpoint = toBreakpointSummary(node);
			String id = idGenerator.generate(breakpoint);
			VariantContextDirectedEvidence result = toVariant(id, node, breakpoint, false); 
			lastHigh = toVariant(id, node, breakpoint, true);
			return result;
		}
		return endOfData();
	}
}
