package au.edu.wehi.idsv;

import java.util.Collection;
import java.util.Iterator;

import au.edu.wehi.idsv.visualisation.TrackedState;
import com.google.common.base.Function;
import com.google.common.collect.AbstractIterator;

import au.edu.wehi.idsv.graph.RectangleGraphMaximalCliqueIterator;
import au.edu.wehi.idsv.graph.RectangleGraphNode;
import au.edu.wehi.idsv.graph.RectangleGraphNodeMergingIterator;
import au.edu.wehi.idsv.graph.ScalingHelper;
import au.edu.wehi.idsv.util.MathUtil;
import au.edu.wehi.idsv.util.WindowedSortingIterator;
import au.edu.wehi.idsv.vcf.VcfInfoAttributes;
import au.edu.wehi.idsv.vcf.VcfSvConstants;
import htsjdk.samtools.util.Log;
/**
 * Maximal clique summary evidence iterator 
 * @author Daniel Cameron
 *
 */
public class MaximalEvidenceCliqueIterator extends AbstractIterator<VariantContextDirectedBreakpoint> implements TrackedState {
	private static final Log log = Log.getInstance(MaximalEvidenceCliqueIterator.class);
	public static final String BREAKPOINT_ID_SUFFIX_HIGH = "h";
	public static final String BREAKPOINT_ID_SUFFIX_LOW = "o";
	private VariantContextDirectedBreakpoint lastHigh = null;
	private final BreakendDirection targetLowDir;
	private final BreakendDirection targetHighDir;
	private final RectangleGraphMaximalCliqueIterator calc;
	private final ProcessingContext context;
	private final VariantIdGenerator idGenerator;
	public MaximalEvidenceCliqueIterator(ProcessingContext processContext, Iterator<DirectedEvidence> evidenceIt, BreakendDirection lowDir, BreakendDirection highDir, VariantIdGenerator idGenerator) {
		this.context = processContext;
		this.calc = new RectangleGraphMaximalCliqueIterator(
						// collapse evidence at the same location to a single node
						new RectangleGraphNodeMergingIterator(RectangleGraphNode.ByStartXYEndXY,
							// make sure nodes to be merged are adjacent in the stream
							new GraphNodeWindowedSortingIterator(context, 1, 
								// convert evidence breakpoints to GraphNodes
								new EvidenceToGraphNodeIterator(evidenceIt))));
		this.targetLowDir = lowDir;
		this.targetHighDir = highDir;
		this.idGenerator = idGenerator;
	}

	private class GraphNodeWindowedSortingIterator extends WindowedSortingIterator<RectangleGraphNode> {
		public GraphNodeWindowedSortingIterator(final GenomicProcessingContext processContext, final int windowSize, final Iterator<RectangleGraphNode> it) {
			super(it, new Function<RectangleGraphNode, Long>() {
				public Long apply(RectangleGraphNode arg) {
					return arg.startX;
				}
			}, windowSize, RectangleGraphNode.ByStartXYEndXY);
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
		if (!(loc instanceof BreakpointSummary)) return null;
		BreakpointSummary bp = (BreakpointSummary)loc;
		if (!bp.isValid(context.getDictionary())) {
			String msg = String.format("Evidence %s has invalid breakpoint %s", e.getEvidenceID(), bp);
			log.error(msg);
			throw new IllegalArgumentException(msg);
		}
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
	private VariantContextDirectedBreakpoint toVariant(String event, RectangleGraphNode node, BreakpointSummary breakpoint, boolean isHighBreakend) {
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(context);
		builder.attribute(VcfSvConstants.BREAKEND_EVENT_ID_KEY, event);
		builder.attribute(VcfSvConstants.PARTNER_BREAKEND_ID_KEY, event + (isHighBreakend ? BREAKPOINT_ID_SUFFIX_LOW : BREAKPOINT_ID_SUFFIX_HIGH));
		builder.id(event + (isHighBreakend ? BREAKPOINT_ID_SUFFIX_HIGH : BREAKPOINT_ID_SUFFIX_LOW));
		if (isHighBreakend) {
			breakpoint = breakpoint.remoteBreakpoint();
		}
		builder.breakpoint(breakpoint, "");
		long scaledWeight = node.weight;
		double weight = ScalingHelper.toUnscaledWeight(scaledWeight);
		builder.phredScore(weight);
		builder.attribute(VcfInfoAttributes.CALLED_QUAL, weight);
		VariantContextDirectedBreakpoint v = (VariantContextDirectedBreakpoint)builder.make();
		assert(v != null);
		return v;
	}
	private BreakpointSummary toBreakpointSummary(RectangleGraphNode node) {
		int start = context.getLinear().getReferencePosition(node.startX);
		int end = context.getLinear().getReferencePosition(node.endX);
		int start2 = context.getLinear().getReferencePosition(node.startY);
		int end2 = context.getLinear().getReferencePosition(node.endY);
		BreakpointSummary breakpoint = new BreakpointSummary(
				context.getLinear().getReferenceIndex(node.startX),
				targetLowDir,
				MathUtil.average(start, end),
				start,
				end,
				context.getLinear().getReferenceIndex(node.startY),
				targetHighDir,
				MathUtil.average(start2, end2),
				start2,
				end2);
		// sanity check that the resultant breakpoint makes sense
		assert(breakpoint.isValid(context.getDictionary()));
		return breakpoint;
	}
	@Override
	protected VariantContextDirectedBreakpoint computeNext() {
		if (lastHigh != null) {
			VariantContextDirectedBreakpoint result = lastHigh;
			lastHigh = null;
			return result;
		}
		if (calc.hasNext()) {
			RectangleGraphNode node = calc.next();
			BreakpointSummary breakpoint = toBreakpointSummary(node);
			String id = idGenerator.generate(breakpoint);
			VariantContextDirectedBreakpoint result = toVariant(id, node, breakpoint, false); 
			lastHigh = toVariant(id, node, breakpoint, true);
			return result;
		}
		return endOfData();
	}

	@Override
	public String[] trackedNames() {
		return calc.trackedNames();
	}

	@Override
	public Object[] trackedState() {
		return calc.trackedState();
	}

	@Override
	public Collection<TrackedState> trackedObjects() {
		return calc.trackedObjects();
	}
}
