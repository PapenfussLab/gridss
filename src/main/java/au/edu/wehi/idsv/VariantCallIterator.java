package au.edu.wehi.idsv;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.Supplier;

import au.edu.wehi.idsv.configuration.VisualisationConfiguration;
import au.edu.wehi.idsv.visualisation.PositionalDeBruijnGraphTracker;
import au.edu.wehi.idsv.visualisation.StateTracker;
import au.edu.wehi.idsv.visualisation.TrackedState;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import htsjdk.samtools.util.Log;
import org.apache.commons.lang3.tuple.Pair;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
/**
 * Calls breakpoints from the given evidence
 * 
 * @author Daniel Cameron
 */
public class VariantCallIterator implements CloseableIterator<VariantContextDirectedEvidence> {
	private static final Log log = Log.getInstance(VariantCallIterator.class);
	private static final List<Pair<BreakendDirection, BreakendDirection>> DIRECTION_ORDER = ImmutableList.of(
			Pair.of(BreakendDirection.Forward, BreakendDirection.Forward),
			Pair.of(BreakendDirection.Forward, BreakendDirection.Backward),
			Pair.of(BreakendDirection.Backward, BreakendDirection.Forward),
			Pair.of(BreakendDirection.Backward, BreakendDirection.Backward),
			Pair.of(BreakendDirection.Forward, null),
			Pair.of(BreakendDirection.Backward, null));
	private final ProcessingContext processContext;
	private final VariantIdGenerator idGenerator;
	private final Supplier<Iterator<DirectedEvidence>> iteratorGenerator;
	private final QueryInterval[] filterInterval;
	private Iterator<? extends VariantContextDirectedEvidence> currentIterator;
	private Iterator<DirectedEvidence> currentUnderlyingIterator;
	private int currentDirectionOrdinal;
	private VariantContextDirectedEvidence lastElement;
	private StateTracker currentTracker;
	private Collection<TrackedState> currentTrackedObjects;
	private PositionTrackedState currentPositionTrackedState = new PositionTrackedState();

	public VariantCallIterator(ProcessingContext processContext, Iterable<DirectedEvidence> evidence) throws InterruptedException {
		this.processContext = processContext;
		this.idGenerator = new SequentialIdGenerator("gridss");
		this.iteratorGenerator = () -> evidence.iterator();
		this.filterInterval = null;
		this.currentDirectionOrdinal = 0;
		reinitialiseIterator();
	}
	public VariantCallIterator(AggregateEvidenceSource source) {
		this.processContext = source.getContext();
		this.idGenerator = new SequentialIdGenerator("gridss");
		this.iteratorGenerator = () -> source.iterator();
		this.filterInterval = null;
		this.currentDirectionOrdinal = 0;
		reinitialiseIterator();
	}
	public VariantCallIterator(AggregateEvidenceSource source, QueryInterval[] interval, int intervalNumber) {
		this.processContext = source.getContext();
		this.idGenerator = new SequentialIdGenerator(String.format("gridss%d_", intervalNumber));
		int expandBy = source.getMaxConcordantFragmentSize() + 1;
		QueryInterval[] expanded = QueryIntervalUtil.padIntervals(processContext.getDictionary(), interval, expandBy);
		this.iteratorGenerator = () -> source.iterator(expanded);
		this.filterInterval = interval;
		this.currentDirectionOrdinal = 0;
		reinitialiseIterator();
	}
	private void reinitialiseIterator() {	
		assert(currentIterator == null || !currentIterator.hasNext());
		CloserUtil.close(currentIterator);
		CloserUtil.close(currentUnderlyingIterator);
		if (currentDirectionOrdinal >= DIRECTION_ORDER.size()) return;
		currentUnderlyingIterator = iteratorGenerator.get();
		Pair<BreakendDirection, BreakendDirection> direction = DIRECTION_ORDER.get(currentDirectionOrdinal);
		if (direction.getRight() != null) {
			currentIterator = new MaximalEvidenceCliqueIterator(
					processContext,
					currentUnderlyingIterator,
					direction.getLeft(),
					direction.getRight(),
					idGenerator);
		} else {
			if (processContext.getVariantCallingParameters().callBreakends) {
				currentIterator = new BreakendMaximalEvidenceCliqueIterator(
						processContext,
						currentUnderlyingIterator,
						direction.getLeft(),
						idGenerator);
			} else {
				if (currentIterator != null) {
					try {
						currentTracker.close();
					} catch (IOException e) {
						log.debug("Telemetry failure", e);
					}
					currentTracker = null;
					currentTrackedObjects = null;
				}
				currentIterator = null;
			}
		}
		if (currentIterator != null && processContext.getConfig().getVisualisation().maxCliqueTelemetry && currentIterator instanceof TrackedState) {
			TrackedState ts = (TrackedState)currentIterator;
			String positionComponent = (filterInterval == null || filterInterval.length == 0) ? "" : String.format("_%s_%d",
					processContext.getDictionary().getSequence(filterInterval[0].referenceIndex).getSequenceName(),
					filterInterval[0].start);
			String filename = String.format("maxclique%s_%s%s.csv", positionComponent, direction.getLeft(), direction.getRight());
			File file = new File(processContext.getConfig().getVisualisation().directory, filename);
			try {
				currentTracker = new StateTracker(file);
				currentTrackedObjects = Lists.newArrayList(Iterables.concat(ts.trackedObjects(), currentPositionTrackedState.trackedObjects()));
				currentTracker.writeHeader(currentTrackedObjects);
			} catch (IOException e) {
				log.debug("Telemetry failure", e);
			}
		}
		if (currentIterator != null && filterInterval != null) {
			currentIterator = Iterators.filter(currentIterator, v -> {
				if (v instanceof DirectedBreakpoint) {
					BreakpointSummary bs = ((DirectedBreakpoint)v).getBreakendSummary();
					return QueryIntervalUtil.overlaps(filterInterval, bs.referenceIndex, bs.start) || 
							QueryIntervalUtil.overlaps(filterInterval, bs.referenceIndex2, bs.start2);
				} else {
					BreakendSummary be = v.getBreakendSummary();
					return QueryIntervalUtil.overlaps(filterInterval, be.referenceIndex, be.start);
				}
			});
		}
	}
	@Override
	public boolean hasNext() {
		if (currentDirectionOrdinal >= DIRECTION_ORDER.size()) return false;
		while (currentIterator != null && currentDirectionOrdinal < DIRECTION_ORDER.size() && !currentIterator.hasNext()) {
			currentDirectionOrdinal++;
			reinitialiseIterator();
		}
		return currentIterator != null && currentIterator.hasNext();
	}
	@Override
	public VariantContextDirectedEvidence next() {
		if (!hasNext()) throw new NoSuchElementException();
		lastElement = currentIterator.next();
		if (currentTracker != null) {
			try {
				currentTracker.track(currentTrackedObjects);
			} catch (IOException e) {
				log.debug("Telemetry failure", e);
			}
		}
		return lastElement;
	}
	@Override
	public void close() {
		CloserUtil.close(currentIterator);
		CloserUtil.close(currentUnderlyingIterator);
	}

	private class PositionTrackedState implements TrackedState {
		@Override
		public String[] trackedNames() {
			return new String[] {
					"chr",
					"pos",
			};
		}

		@Override
		public Object[] trackedState() {
			return new Object[] {
					lastElement == null ? "" : lastElement.getContig(),
					lastElement == null ? 0 : lastElement.getStart(),
			};
		}

		@Override
		public Collection<TrackedState> trackedObjects() {
			return ImmutableList.of(this);
		}
	}
}
 