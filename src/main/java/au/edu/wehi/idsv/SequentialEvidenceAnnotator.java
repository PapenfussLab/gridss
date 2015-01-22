package au.edu.wehi.idsv;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;
import java.util.PriorityQueue;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Iterators;
import com.google.common.collect.Ordering;
import com.google.common.collect.PeekingIterator;
import com.google.common.collect.TreeMultimap;

/**
 * Annotates sorted breakends with the evidence supporting the call
 * 
 * @author Daniel Cameron
 *
 */
public class SequentialEvidenceAnnotator extends AbstractIterator<VariantContextDirectedEvidence> implements BreakendAnnotator {
	private final ProcessingContext context;
	private final int maxCallRange;
	private final boolean assignEvidenceToSingleBreakpoint;
	private final Iterator<DirectedEvidence> evidenceIt;
	private final PeekingIterator<VariantContextDirectedEvidence> callIt;
	private final TreeMultimap<Long, ActiveVariant> activeVariantStart = TreeMultimap.<Long, ActiveVariant>create(Ordering.natural(), ActiveVariantByEndStart);
	private final PriorityQueue<VariantContextDirectedEvidence> outputBuffer = new PriorityQueue<VariantContextDirectedEvidence>(1024, VariantContextDirectedEvidence.ByLocationStart);
	private final EvidenceToCsv dump;
	/**
	 * Incorporates the next evidence record in the called variants
	 */
	private void processNextEvidence() {
		boolean evidenceCalled = false;
		DirectedEvidence evidence = evidenceIt.next();
		BreakendSummary bs = evidence.getBreakendSummary();
		bs = context.getVariantCallingParameters().withMargin(context, bs);
		long evidenceStart = context.getLinear().getStartLinearCoordinate(bs);
		long evidenceEnd = context.getLinear().getEndLinearCoordinate(bs);
		flushCallsStartBefore(evidenceStart - maxCallRange - 1);
		addCallsStartBefore(evidenceEnd + 1);
		if (assignEvidenceToSingleBreakpoint) {
			float bestScore = Float.MIN_VALUE;
			ActiveVariant best = null;
			for (Collection<ActiveVariant> cv : activeVariantStart.asMap().subMap(evidenceStart, true, evidenceEnd, true).values()) {
				for (ActiveVariant v : cv) {
					if (v.location.overlaps(bs)) {
						if (v.score > bestScore) {
							bestScore = v.score;
							best = v;
						}
					}
				}
			}
			if (best != null) {
				best.attributeEvidence(evidence);
				evidenceCalled = true;
			}
		} else {
			for (Collection<ActiveVariant> cv : activeVariantStart.asMap().subMap(evidenceStart, true, evidenceEnd, true).values()) {
				for (ActiveVariant v : cv) {
					if (v.location.overlaps(bs)) {
						v.attributeEvidence(evidence);
						evidenceCalled = true;
					}
				}
			}
		}
		if (!evidenceCalled && dump != null) {
			// evidence stranded and does not provide support for any call
			// write out now before we drop it
			dump.writeEvidence(evidence, null);
		}
	}
	/**
	 * Adds variant calls starting before the given position to the active set
	 * @param pos linear genomic position
	 */
	private void addCallsStartBefore(long pos) {
		while (callIt.hasNext() && context.getLinear().getStartLinearCoordinate(callIt.peek().getBreakendSummary()) < pos) {
			ActiveVariant av = new ActiveVariant(callIt.next());
			activeVariantStart.put(av.startLocation, av);
		}
	}
	/**
	 * Flushes all variant calls starting before the given linear genomic position 
	 * @param pos linear genomic position
	 */
	private void flushCallsStartBefore(long pos) {
		Iterator<Entry<Long, ActiveVariant>> it = activeVariantStart.entries().iterator();
		while (it.hasNext()) {
			Entry<Long, ActiveVariant> entry = it.next();
			ActiveVariant av = entry.getValue();
			if (av.startLocation >= pos) {
				return;
			}
			it.remove();
			VariantContextDirectedEvidence result = av.callVariant();
			outputBuffer.add(result);
		}
	}
	private class ActiveVariant {
		public final long startLocation;
		public final long endLocation;
		public final BreakendSummary location;
		public final float score;
		private final StructuralVariationCallBuilder builder;
		private List<DirectedEvidence> evidence = new ArrayList<DirectedEvidence>();
		public ActiveVariant(VariantContextDirectedEvidence call) {
			this.builder = new StructuralVariationCallBuilder(context, call);
			this.score = (float)call.getPhredScaledQual();
			this.location = call.getBreakendSummary();
			this.startLocation = context.getLinear().getStartLinearCoordinate(this.location);
			this.endLocation = context.getLinear().getEndLinearCoordinate(this.location);
		}
		public void attributeEvidence(DirectedEvidence e) {
			evidence.add(e);
			builder.addEvidence(e);
		}
		public VariantContextDirectedEvidence callVariant() {
			VariantContextDirectedEvidence call = builder.make();
			if (dump != null) {
				for (DirectedEvidence e : evidence) {
					dump.writeEvidence(e, call);
				}
			}
			return builder.make();
		}
	}
	private static final Ordering<ActiveVariant> ActiveVariantByEndStart = new Ordering<ActiveVariant>() {
		@Override
		public int compare(ActiveVariant arg0, ActiveVariant arg1) {
			 return ComparisonChain.start()
				.compare(arg0.endLocation, arg1.endLocation)
		        .compare(arg0.startLocation, arg1.startLocation)
		        .result();
		}
	};
	public SequentialEvidenceAnnotator(
			ProcessingContext context,
			Iterator<VariantContextDirectedEvidence> calls,
			Iterator<DirectedEvidence> evidence,
			int maxCallWindowSize,
			boolean assignEvidenceToSingleBreakpoint,
			EvidenceToCsv dump) {
		this.context = context;
		this.maxCallRange = maxCallWindowSize;
		this.callIt = Iterators.peekingIterator(calls);
		this.evidenceIt = evidence;
		this.assignEvidenceToSingleBreakpoint = assignEvidenceToSingleBreakpoint;
		this.dump = dump;
	}
	@Override
	protected VariantContextDirectedEvidence computeNext() {
		while (outputBuffer.isEmpty() && evidenceIt.hasNext()) {
			processNextEvidence();
		}
		if (!evidenceIt.hasNext()) {
			addCallsStartBefore(Long.MAX_VALUE);
			flushCallsStartBefore(Long.MAX_VALUE);
		}
		if (!outputBuffer.isEmpty()) {
			return outputBuffer.poll();
		}
		return endOfData();
	}
}
