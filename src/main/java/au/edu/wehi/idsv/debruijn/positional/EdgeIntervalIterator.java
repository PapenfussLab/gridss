package au.edu.wehi.idsv.debruijn.positional;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.PriorityQueue;

/**
 * Iterates over all edge sets for a KmerPathNode
 * 
 *  Note: to eliminate costly collection copy operations, Iterator<T> is not implemented.
 *  Only the current element can be accessed. Collection properties return  
 *  
 * @author Daniel Cameron
 *
 */
public class EdgeIntervalIterator {
	public Collection<KmerPathNode> currentNext() { return validNext; }
	public Collection<KmerPathNode> currentPrev() { return validPrev; }
	/**
	 * First kmer start position
	 * @return
	 */
	public int currentFirstKmerStartPosition() { return start; }
	/**
	 * First kmer end position
	 * @return
	 */
	public int currentFirstKmerEndPosition() { return end; }
	
	private final PriorityQueue<KmerPathNode> validNext = new PriorityQueue<KmerPathNode>(16, KmerPathNode.ByFirstKmerEndPosition);
	private final PriorityQueue<KmerPathNode> validPrev = new PriorityQueue<KmerPathNode>(16, KmerPathNode.ByEndPosition);
	private final int maxEnd;
	private final int length;
	private final List<KmerPathNode> next;
	private final List<KmerPathNode> prev;
	private int nextOffset = 0;
	private int prevOffset = 0;
	private int start = Integer.MIN_VALUE;
	private int end = Integer.MIN_VALUE;
	public EdgeIntervalIterator(KmerPathNode node) {
		this(node.startPosition(), node.endPosition(), node.length(), node.next(), node.prev(), false);
	}
	public EdgeIntervalIterator(int start, int end, int length, List<KmerPathNode> next, List<KmerPathNode> prev, boolean preSorted) {
		if (!preSorted) {
			Collections.sort(next, KmerPathNode.ByFirstKmerStartPosition);
			Collections.sort(prev, KmerPathNode.ByStartPosition);
		}
		this.next = next;
		this.prev = prev;
		this.maxEnd = end;
		this.length = length;
		advanceTo(start);
	}
	private void advanceTo(int startPos) {
		this.start = startPos;
		while (nextOffset < next.size() && next.get(nextOffset).startPosition(0) <= start + length) {
			validNext.add(next.get(nextOffset++));
		}
		while (prevOffset < prev.size() && prev.get(prevOffset).startPosition() <= start - 1) {
			validNext.add(prev.get(prevOffset++));
		}
		this.end = maxEnd;
		if (!validNext.isEmpty()) {
			this.end = Math.min(maxEnd, validNext.peek().endPosition(0) - length);
		}
		if (!validPrev.isEmpty()) {
			this.end = Math.min(maxEnd, validPrev.peek().endPosition() + 1);
		}
	}
	public void advance() {
		if (end >= maxEnd) {
			start = Integer.MAX_VALUE;
			end = Integer.MAX_VALUE;
		} else {
			// next position is the first change encountered
			int nextStartPos = Integer.MAX_VALUE;
			if (!validNext.isEmpty()) {
				nextStartPos = Math.min(nextStartPos, validNext.peek().endPosition(0) - length + 1);
			}
			if (!validPrev.isEmpty()) {
				nextStartPos = Math.min(nextStartPos, validPrev.peek().endPosition() + 1 + 1);
			}
			if (nextOffset < next.size() && next.get(nextOffset).startPosition(0) <= start + length) {
				nextStartPos = Math.min(nextStartPos, next.get(nextOffset).startPosition(0) - length);
			}
			if (prevOffset < prev.size() && prev.get(prevOffset).startPosition() <= start - 1) {
				nextStartPos = Math.min(nextStartPos, prev.get(prevOffset).startPosition() + 1);
			}
			// flush out neighbours that end
			while (!validNext.isEmpty() && nextStartPos == validNext.peek().endPosition(0) - length + 1) {
				validNext.poll();
			}
			while (!validPrev.isEmpty() && nextStartPos == validPrev.peek().endPosition() + 1 + 1) {
				validPrev.poll();
			}
			advanceTo(nextStartPos);
		}
	}
	public boolean isValid() {
		return end - start > 0;
	}
}