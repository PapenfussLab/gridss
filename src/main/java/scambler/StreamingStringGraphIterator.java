package scambler;

import java.util.ArrayDeque;
import java.util.Deque;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.PriorityQueue;

import au.edu.wehi.idsv.LinearGenomicCoordinate;
import htsjdk.samtools.SAMRecord;

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
		toTransform.push(r);
		currentPosition = lgc.getStartLinearCoordinate(record);
		// transform records
		while (getWindowPosition(toTransform.peek()) < currentTransformPosition()) {
			transform(toTransform.pop());
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
				transform(toTransform.pop());
			}
		}
	}
	@Override
	public boolean hasNext() {
		ensureNext();
		return !output.isEmpty();
	}
	@Override
	public SgNode next() {
		ensureNext();
		if (output.isEmpty()) {
			throw new NoSuchElementException();
		}
		SgNode node = output.poll();
		if (node.read.getEndNode() == node) {
			// Each read creates both a start and end node.
			// Remove the read from the lookup when we emit the
			// end node
			lookup.remove(node.read);
		}
		return node;
	}
}