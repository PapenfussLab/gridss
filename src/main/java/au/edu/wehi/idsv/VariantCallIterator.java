package au.edu.wehi.idsv;

import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.function.Supplier;

import org.apache.commons.lang3.tuple.Pair;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.util.CloserUtil;
/**
 * Calls breakpoints from the given evidence
 * 
 * @author Daniel Cameron
 */
public class VariantCallIterator implements Iterator<VariantContextDirectedEvidence> {
	private static final List<Pair<BreakendDirection, BreakendDirection>> DIRECTION_ORDER = ImmutableList.of(
			Pair.of(BreakendDirection.Forward, BreakendDirection.Forward),
			Pair.of(BreakendDirection.Forward, BreakendDirection.Backward),
			Pair.of(BreakendDirection.Backward, BreakendDirection.Forward),
			Pair.of( BreakendDirection.Backward, BreakendDirection.Backward));
	private final ProcessingContext processContext;
	private final VariantIdGenerator idGenerator;
	private final Supplier<Iterator<DirectedEvidence>> iteratorGenerator;
	private final QueryInterval filterInterval;
	private Iterator<VariantContextDirectedEvidence> currentIterator;
	private int currentDirectionOrdinal;
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
	public VariantCallIterator(AggregateEvidenceSource source, QueryInterval interval, int intervalNumber) {
		int expandBy = source.getMaxConcordantFragmentSize() + 1;
		QueryInterval inputInterval = new QueryInterval(interval.referenceIndex, interval.start - expandBy, interval.end + expandBy);
		this.processContext = source.getContext();
		this.idGenerator = new SequentialIdGenerator(String.format("gridss%d_", intervalNumber));
		this.iteratorGenerator = () -> source.iterator(inputInterval);
		this.filterInterval = interval;
		this.currentDirectionOrdinal = 0;
		reinitialiseIterator();
	}
	private void reinitialiseIterator() {
		assert(currentIterator == null || !currentIterator.hasNext());
		CloserUtil.close(currentIterator);
		if (currentDirectionOrdinal >= DIRECTION_ORDER.size()) return; 
		currentIterator = new MaximalEvidenceCliqueIterator(
				processContext,
				iteratorGenerator.get(),
				DIRECTION_ORDER.get(currentDirectionOrdinal).getLeft(),
				DIRECTION_ORDER.get(currentDirectionOrdinal).getRight(),
				idGenerator);
		if (filterInterval != null) {
			currentIterator = Iterators.filter(currentIterator, v -> {
				BreakendSummary bs = v.getBreakendSummary();
				return bs.start >= filterInterval.start && bs.start <= filterInterval.end; 
			});
		}
	}
	@Override
	public boolean hasNext() {
		if (currentDirectionOrdinal >= DIRECTION_ORDER.size()) return false;
		if (!currentIterator.hasNext()) {
			currentDirectionOrdinal++;
			reinitialiseIterator();
			return hasNext();
		}
		return true;
	}
	@Override
	public VariantContextDirectedEvidence next() {
		if (!hasNext()) throw new NoSuchElementException();
		return currentIterator.next();
	}
}
 