package scambler;

import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.LinearGenomicCoordinate;
import htsjdk.samtools.SAMRecord;

import java.util.*;

public class StreamingStringGraphIterator implements Iterator<SgNode> {
	private final OverlapLookup lookup;
	private final Deque<Read> toTransform = new ArrayDeque<>();
	private final long windowSize;
	private final int maxAssemblyOverlapDistance;
	private final Iterator<SAMRecord> it;
	private final LinearGenomicCoordinate lgc;
	private final PriorityQueue<SgNode> output = new PriorityQueue<>(SgNode.ByInferredPosition);
	private long currentPosition;
	public StreamingStringGraphIterator(int minReadOverlap, int maxAssemblyOverlapDistance, int maxReadLength, Iterator<SAMRecord> it, LinearGenomicCoordinate lgc) {
		this.lookup = new OverlapLookup(minReadOverlap);
		// Need maxReadLength since soft clipping a read moves the inferred position to before the sort order position
		this.windowSize = maxAssemblyOverlapDistance + maxReadLength;
		this.maxAssemblyOverlapDistance = maxAssemblyOverlapDistance;
		this.lgc = lgc;
		this.it = it;
	}
	private void load(SAMRecord record) {
		Read r = Read.create(lgc, record);
		output.add(r.getStartNode());
		output.add(r.getEndNode());
		lookup.add(r);
		toTransform.addLast(r);
		long nextPosition = lgc.getStartLinearCoordinate(record);
		assert(nextPosition >= currentPosition);
		currentPosition = nextPosition;
		// transform records
		while (getWindowPosition(toTransform.peekFirst()) < currentTransformPosition()) {
			transform(toTransform.removeFirst());
		}
	}
	private void transform(Read r) {
		assert(currentPosition > getWindowPosition(r) + windowSize);
		// calculate overlaps
		for (Overlap o : lookup.successors(r)) {
			// only consider overlaps which are within the assembly window
			if (o.deviationFromAlignment() <= maxAssemblyOverlapDistance) {
				SgEdge.create(o);
			}
		}
	}
	private long getWindowPosition(Read r) {
		// TODO: update to inferred position interval
		return lgc.getLinearCoordinate(r.getRead().getReferenceIndex(), (r.getRead().getUnclippedEnd() + r.getRead().getUnclippedStart()) / 2);
	}
	private long currentTransformPosition() { return currentPosition - windowSize; }
	private long currentEmitPosition() { return currentPosition - 2 * windowSize; }
	private void ensureNext() {
		while (it.hasNext() && (output.isEmpty() || getWindowPosition(output.peek().read) > currentEmitPosition())) {
			SAMRecord r = it.next();
			load(r);
		}
		if (!it.hasNext()) {
			currentPosition = Long.MAX_VALUE;
			while (!toTransform.isEmpty()) {
				transform(toTransform.removeFirst());
			}
		}
	}
	private void ensureNextNonOverlapping() {
		ensureNext();
		while (!output.isEmpty() && isNonoverlapping(output.peek())) {
			output.poll();
			ensureNext();
		}
	}
	private boolean isNonoverlapping(SgNode node) {
		return node.in.size() == 0 && node.out.size() == 0;
	}
	@Override
	public boolean hasNext() {
		ensureNextNonOverlapping();
		return !output.isEmpty();
	}
	@Override
	public SgNode next() {
		ensureNextNonOverlapping();
		if (output.isEmpty()) {
			throw new NoSuchElementException();
		}
		SgNode node = output.poll();
		assert(toTransform.isEmpty() || getWindowPosition(toTransform.peekFirst()) > node.inferredPosition + windowSize);
		if (node.read.getEndNode() == node) {
			// Each read creates both a start and end node.
			// Remove the read from the lookup when we emit the
			// end node
			lookup.remove(node.read);
			if (Defaults.SANITY_CHECK_ASSEMBLY_GRAPH) {
				// make sure we never emit a node that could potentially overlap with
				assert(toTransform.stream().allMatch(r -> getWindowPosition(r) > node.inferredPosition + windowSize));
			}
		}
		return node;
	}
}