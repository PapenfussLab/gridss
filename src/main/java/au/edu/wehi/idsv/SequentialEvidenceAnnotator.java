package au.edu.wehi.idsv;

import java.util.Iterator;
import java.util.List;
import java.util.concurrent.Executor;

import au.edu.wehi.idsv.util.ParallelTransformIterator;
import au.edu.wehi.idsv.visualisation.TrackedBuffer;

/**
 * Annotates each variant call based on supporting evidence. Both the variant calls and
 * evidence are required to be in-order.
 * @author Daniel Cameron
 *
 */
public class SequentialEvidenceAnnotator extends ParallelTransformIterator<StructuralVariationCallBuilder, VariantContextDirectedEvidence> implements TrackedBuffer {
	/**
	 * Creates an ordered evidence annotation
	 * @param context processing context
	 * @param calls variant calls ordered by position
	 * @param evidence evidence order by breakend position
	 * @param maxCallWindowSize
	 * @param assignEvidenceToSingleBreakpoint uniquely assign evidence to only the highest scoring call
	 * @param lookahead number of records to annotate in parallel.
	 * @param threadpool thread pool used to perform annotation
	 */
	public SequentialEvidenceAnnotator(
			ProcessingContext context,
			Iterator<? extends VariantContextDirectedEvidence> calls,
			Iterator<? extends DirectedEvidence> evidence,
			int maxCallWindowSize,
			boolean assignEvidenceToSingleBreakpoint,
			int lookahead,
			Executor threadpool) {
		super(new SequentialEvidenceAllocator(context, calls, evidence, maxCallWindowSize, assignEvidenceToSingleBreakpoint),
			call -> call.make(), lookahead, threadpool);
	}

	@Override
	public void setTrackedBufferContext(String context) {
		((TrackedBuffer)it).setTrackedBufferContext(context);
	}

	@Override
	public List<NamedTrackedBuffer> currentTrackedBufferSizes() {
		return ((TrackedBuffer)it).currentTrackedBufferSizes();
	}
}
