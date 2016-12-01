package au.edu.wehi.idsv;

import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.concurrent.Executor;
import java.util.stream.Collectors;

import au.edu.wehi.idsv.SequentialEvidenceAllocator.VariantEvidenceSupport;
import au.edu.wehi.idsv.util.ParallelTransformIterator;
import au.edu.wehi.idsv.vcf.VcfAttributes;
import au.edu.wehi.idsv.visualisation.TrackedBuffer;
import htsjdk.samtools.util.Log;

/**
 * Annotates each variant call based on supporting evidence. Both the variant calls and
 * evidence are required to be in-order.
 * @author Daniel Cameron
 *
 */
public class SequentialEvidenceAnnotator extends ParallelTransformIterator<VariantEvidenceSupport, VariantContextDirectedEvidence> implements TrackedBuffer {
	private static final Log log = Log.getInstance(SequentialEvidenceAnnotator.class);
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
			call -> make(context, call), lookahead, threadpool);
	}
	private static VariantContextDirectedEvidence make(ProcessingContext context, VariantEvidenceSupport ves) {
		try {
			StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(context, ves.variant);
			Object eid = ves.variant.getAttribute(VcfAttributes.EVIDENCE_ID.attribute());
			if (eid == null) {
				for (DirectedEvidence e : ves.support) {
					builder.addEvidence(e);
				}
			} else {
				// only add evidence that has been assigned
				Set<String> allocatedEvidence = ves.variant.getAttributeAsList(VcfAttributes.EVIDENCE_ID.attribute())
						.stream()
						.map(s -> (String)s)
						.collect(Collectors.toSet());
				for (DirectedEvidence e : ves.support) {
					if (allocatedEvidence.contains(e.getEvidenceID())) {
						builder.addEvidence(e);
					}
				}
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
