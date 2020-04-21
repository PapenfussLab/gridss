package au.edu.wehi.idsv;

import au.edu.wehi.idsv.SequentialEvidenceAllocator.VariantEvidenceSupport;
import au.edu.wehi.idsv.util.ParallelTransformIterator;
import au.edu.wehi.idsv.visualisation.TrackedBuffer;
import htsjdk.samtools.util.Log;

import java.util.Iterator;
import java.util.List;
import java.util.concurrent.Executor;

/**
 * Annotates each variant call based on supporting evidence. Both the variant calls and
 * evidence are required to be in-order.
 * @author Daniel Cameron
 *
 */
public class SequentialEvidenceAnnotator extends ParallelTransformIterator<VariantEvidenceSupport, VariantContextDirectedEvidence> implements TrackedBuffer {
	private static final Log log = Log.getInstance(SequentialEvidenceAnnotator.class);
	private final ProcessingContext context;
	private final CalledBreakpointPositionLookup lookup = new CalledBreakpointPositionLookup();

	/**
	 * Creates an ordered evidence annotation
	 * @param context processing context
	 * @param calls variant calls ordered by position
	 * @param reads reads order by breakend position
	 * @param assemblies assemblies order by breakend position
	 * @param maxCallWindowSize
	 * @param assignEvidenceToSingleBreakpoint uniquely assign evidence to only the highest scoring call
	 * @param lookahead number of records to annotate in parallel.
	 * @param threadpool thread pool used to perform annotation
	 */
	public SequentialEvidenceAnnotator(
			ProcessingContext context,
			Iterator<? extends VariantContextDirectedEvidence> calls,
			Iterator<? extends DirectedEvidence> reads,
			Iterator<? extends DirectedEvidence> assemblies,
			int maxCallWindowSize,
			boolean assignEvidenceToSingleBreakpoint,
			int lookahead,
			Executor threadpool) {
		super(new SequentialEvidenceAllocator(context, calls, reads, assemblies, maxCallWindowSize, assignEvidenceToSingleBreakpoint),
				null, lookahead, threadpool);
		this.context = context;
	}

	@Override
	public VariantContextDirectedEvidence transform(VariantEvidenceSupport ves) {
		try {
			StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(context, lookup, ves.variant);
			for (DirectedEvidence e : ves.support) {
				builder.addEvidence(e);
			}
			VariantContextDirectedEvidence evidence = builder.make();
			return evidence;
		} catch (Exception e) {
			log.error("Error annotating ", ves.variant.getID());
			throw e;
		}
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
