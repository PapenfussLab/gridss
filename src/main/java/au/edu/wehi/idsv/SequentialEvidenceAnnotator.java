package au.edu.wehi.idsv;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

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
	private final PeekingIterator<? extends DirectedEvidence> evidenceIt;
	private final Iterator<? extends VariantContextDirectedEvidence> callIt;
	private final ArrayDeque<ActiveVariant> variantBuffer = new ArrayDeque<ActiveVariant>();
	private final EvidenceToCsv dump;
	private class ActiveVariant {
		public final long startLocation;
		//public final long endLocation;
		public final BreakendSummary location;
		public final float score;
		private final StructuralVariationCallBuilder builder;
		private List<DirectedEvidence> evidence = new ArrayList<DirectedEvidence>();
		public ActiveVariant(VariantContextDirectedEvidence call) {
			this.builder = new StructuralVariationCallBuilder(context, call);
			this.score = (float)call.getPhredScaledQual();
			assert(this.score >= 0); // variant must have score set
			this.location = call.getBreakendSummary();
			this.startLocation = context.getLinear().getStartLinearCoordinate(this.location);
			//this.endLocation = context.getLinear().getEndLinearCoordinate(this.location);
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
	public SequentialEvidenceAnnotator(
			ProcessingContext context,
			Iterator<? extends VariantContextDirectedEvidence> calls,
			Iterator<? extends DirectedEvidence> evidence,
			int maxCallWindowSize,
			boolean assignEvidenceToSingleBreakpoint,
			EvidenceToCsv dump) {
		this.context = context;
		this.maxCallRange = maxCallWindowSize;
		this.callIt = calls;
		this.evidenceIt = Iterators.peekingIterator(evidence);
		this.assignEvidenceToSingleBreakpoint = assignEvidenceToSingleBreakpoint;
		this.dump = dump;
	}
	private void buffer(VariantContextDirectedEvidence variant) {
		variantBuffer.add(new ActiveVariant(variant));
	}
	@Override
	protected VariantContextDirectedEvidence computeNext() {
		if (variantBuffer.isEmpty()) {
			if (!callIt.hasNext()) {
				if (Defaults.PERFORM_ITERATOR_SANITY_CHECKS) {
					// advance so we can check it was correct, even though
					// we have no more calls to make so this doesn't
					// actually need to be done
					while (evidenceIt.hasNext()) evidenceIt.next();
				}
				return endOfData();
			}
			buffer(callIt.next());
		}
		ActiveVariant variant = variantBuffer.peek();
		bufferVariantsBefore(variant.startLocation + 2 * (maxCallRange + 1));
		processEvidenceBefore(variant.startLocation + maxCallRange + 1);
		variantBuffer.poll();
		return variant.callVariant();
	}
	private void processEvidenceBefore(long position) {
		while (evidenceIt.hasNext() && context.getLinear().getStartLinearCoordinate(evidenceIt.peek().getBreakendSummary()) - context.getVariantCallingParameters().breakendMargin <= position) {
			assignEvidence(evidenceIt.next());
		}
	}
	/**
	 * Assigns the given evidence to the appropriate call
	 * 
	 */
	private void assignEvidence(DirectedEvidence evidence) {
		// TODO: replace overlap implementation with interval tree
		boolean evidenceCalled = false;
		BreakendSummary bs = evidence.getBreakendSummary();
		bs = context.getVariantCallingParameters().withMargin(context, bs);
		long endLocation = context.getLinear().getEndLinearCoordinate(bs);
		if (assignEvidenceToSingleBreakpoint) {
			float bestScore = Float.MIN_VALUE;
			ActiveVariant best = null;
			for (ActiveVariant v : variantBuffer) {
				if (v.startLocation > endLocation) break;
				if (v.location.overlaps(bs)) {
					if (v.score > bestScore) {
						bestScore = v.score;
						best = v;
					}
				}
			}
			if (best != null) {
				best.attributeEvidence(evidence);
				evidenceCalled = true;
			}
		} else {
			for (ActiveVariant v : variantBuffer) {
				if (v.startLocation > endLocation) break;
				if (v.location.overlaps(bs)) {
					v.attributeEvidence(evidence);
					evidenceCalled = true;
				}
			}
		}
		if (!evidenceCalled && dump != null) {
			// evidence does not provide support for any call
			// write out now before we drop it
			dump.writeEvidence(evidence, null);
		}
	}
	private void bufferVariantsBefore(long position) {
		while (callIt.hasNext() && (variantBuffer.isEmpty() || variantBuffer.peekLast().startLocation <= position)) {
			buffer(callIt.next());
		}
	}
}
