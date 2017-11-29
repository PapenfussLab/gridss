package scambler;

import java.util.ArrayDeque;
import java.util.Deque;
import java.util.Iterator;
import java.util.List;

import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.LinearGenomicCoordinate;
import htsjdk.samtools.SAMRecord;
import it.uniroma1.dis.wsngroup.gexf4j.core.Edge;
import it.uniroma1.dis.wsngroup.gexf4j.core.EdgeType;
import it.uniroma1.dis.wsngroup.gexf4j.core.Node;

public class StreamingStringGraphIterator implements Iterator<SgNode> {
	private final OverlapLookup lookup;
	private final Deque<Read> reads = new ArrayDeque<>();
	private final long maxAssemblyOverlapDistance;
	private final Iterator<SAMRecord> it;
	private final LinearGenomicCoordinate lgc;
	private long currentPosition;
	public StreamingStringGraphIterator(int minReadOverlap, int maxAssemblyOverlapDistance, Iterator<SAMRecord> it, LinearGenomicCoordinate lgc) {
		this.lookup = new OverlapLookup(minReadOverlap);
		this.maxAssemblyOverlapDistance = maxAssemblyOverlapDistance;
		this.lgc = lgc;
		this.it = it;
	}
	private void load(SAMRecord record) {
		Read r = Read.create(record);
		reads.push(r);
		lookup.add(r);
	}
	private void transform(Read r) {
		assert(currentPosition > getWindowEndPosition(r) + maxAssemblyOverlapDistance);
		// calculate overlaps
		for (Overlap o : lookup.successors(r)) {
			List<SgEdge> edges = SgEdge.create(o);
		}
	}
	private long getWindowStartPosition(Read r) {
		return lgc.getLinearCoordinate(r.getRead().getReferenceIndex(), r.getRead().getUnclippedStart());
	}
	private long getWindowEndPosition(Read r) {
		return lgc.getLinearCoordinate(r.getRead().getReferenceIndex(), r.getRead().getUnclippedEnd());
	}
	private void unloadBefore(long position) {
		while (!reads.isEmpty()) {
			if (getWindowEndPosition(reads.peek()) < position) {
				Read r = reads.pop();
				lookup.remove(r);
			}
		}
	}
}