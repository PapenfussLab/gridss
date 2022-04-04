package au.edu.wehi.idsv;

import au.edu.wehi.idsv.graph.ScalingHelper;
import au.edu.wehi.idsv.vcf.VcfInfoAttributes;
import au.edu.wehi.idsv.vcf.VcfSvConstants;
import au.edu.wehi.idsv.visualisation.TrackedState;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

import java.util.Collection;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.PriorityQueue;

public class BreakendMaximalEvidenceCliqueIterator implements Iterator<VariantContextDirectedEvidence>, TrackedState {
	public static final String BREAKEND_ID_SUFFIX = "b";
	private final BreakendDirection direction;
	private final ProcessingContext context;
	private final VariantIdGenerator idGenerator;
	private PeekingIterator<DirectedEvidence> it;
	private long activeScore = 0;
	private PriorityQueue<DirectedEvidence> activeByEnd = new PriorityQueue<>(DirectedEvidence.ByEndStart);
	public BreakendMaximalEvidenceCliqueIterator(
			ProcessingContext processContext,
			Iterator<DirectedEvidence> it,
			BreakendDirection direction,
			VariantIdGenerator idGenerator) {
		this.context = processContext;
		this.idGenerator = idGenerator;
		this.direction = direction;
		this.it = Iterators.peekingIterator(
				Iterators.filter(it,
				de -> de.getBreakendSummary().direction == direction &&
				!(de instanceof DirectedBreakpoint) &&
				ScalingHelper.toScaledWeight(de.getBreakendQual()) > 0));
	}
	@Override
	public boolean hasNext() {
		return it.hasNext();
	}
	@Override
	public VariantContextDirectedEvidence next() {
		if (!it.hasNext()) throw new NoSuchElementException();
		LinearGenomicCoordinate lgc = context.getLinear();
		long activeStart = lgc.getStartLinearCoordinate(it.peek().getBreakendSummary());
		// remove evidence whose interval finishes before we start
		while (!activeByEnd.isEmpty() && lgc.getEndLinearCoordinate(activeByEnd.peek().getBreakendSummary()) < activeStart) {
			DirectedEvidence out = activeByEnd.poll();
			activeScore -= ScalingHelper.toScaledWeight(out.getBreakendQual());
		}
		while (it.hasNext() &&
				(activeByEnd.isEmpty() || 
				lgc.getStartLinearCoordinate(it.peek().getBreakendSummary()) <=
				lgc.getEndLinearCoordinate(activeByEnd.peek().getBreakendSummary()))) {
			// this record can be added to our active clique without any removal
			DirectedEvidence de = it.next();
			BreakendSummary bs = de.getBreakendSummary();
			assert(bs.direction == direction);
			float weight = de.getBreakendQual();
			long scaledWeight = ScalingHelper.toScaledWeight(weight);
			assert(scaledWeight > 0);
			activeStart = lgc.getStartLinearCoordinate(bs);
			activeScore += scaledWeight;
			activeByEnd.add(de);
		}
		long activeEnd = lgc.getEndLinearCoordinate(activeByEnd.peek().getBreakendSummary());
		int referenceIndex = lgc.getReferenceIndex(activeStart);
		assert(lgc.getReferenceIndex(activeEnd) == referenceIndex);
		int start = lgc.getReferencePosition(activeStart);
		int end = lgc.getReferencePosition(activeEnd);
		return createRecord(referenceIndex, start, end, ScalingHelper.toUnscaledWeight(activeScore));
	}
	private VariantContextDirectedEvidence createRecord(int referenceIndex, int start, int end, double qual) {
		BreakendSummary breakend = new BreakendSummary(referenceIndex, direction, (start + end) / 2, start, end);
		String id = idGenerator.generate(breakend);
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(context);
		builder.attribute(VcfSvConstants.EVENT_ID_KEY, id);
		builder.id(id + BREAKEND_ID_SUFFIX);
		builder.breakend(breakend, "");
		builder.phredScore(qual);
		builder.attribute(VcfInfoAttributes.CALLED_QUAL, qual);
		VariantContextDirectedEvidence v = (VariantContextDirectedEvidence)builder.make();
		assert(v != null);
		return v;
	}

	@Override
	public String[] trackedNames() {
		return new String[] {"activeByEndSize"};
	}

	@Override
	public Object[] trackedState() {
		return new Object[] {
				activeByEnd == null ? 0 : activeByEnd.size(),
		};
	}

	@Override
	public Collection<TrackedState> trackedObjects() {
		return ImmutableList.of(this);
	}
}
