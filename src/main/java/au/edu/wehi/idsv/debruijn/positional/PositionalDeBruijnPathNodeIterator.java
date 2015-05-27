package au.edu.wehi.idsv.debruijn.positional;

import java.util.Iterator;
import java.util.PriorityQueue;
import java.util.Set;

import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

/**
 * Transforms a start position sorted sequence of KmerNodes to a
 * start position sorted sequence of KmerAggregateNodes
 * @author Daniel Cameron
 *
 */
public class PositionalDeBruijnPathNodeIterator implements Iterator<KmerPathNode> {
	private final PeekingIterator<KmerAggregateNode> underlying;
	private final int maxSupportWidth;
	private PriorityQueue<KmerPathNode> buffer = new PriorityQueue<KmerPathNode>(1024, KmerPathNode.ByFirstKmerStartPosition);
	/**
	 * A PathNode is considered active if the [start, end] interval of the final kmer overlaps
	 * the start position of the last node encountered.
	 *
	 * Active front can be divided into maxSupportWidth sections
	 * - the maxSupportWidth nearest the most recently added KmerAggregateNodes
	 * can have prev or next nodes added
	 * - the maxSupportWidth before that have all adjacent nodes already added to the active set
	 * 
	 * PathNodes are emitted [-2 * maxSupportWidth, -maxSupportWidth] before the underlying iterator position
	 * 
	 */
	private Set<KmerAggregateNode> active;
	public PositionalDeBruijnPathNodeIterator(Iterator<KmerAggregateNode> it, int maxSupportWidth) {
		this.underlying = Iterators.peekingIterator(it);
		this.maxSupportWidth = maxSupportWidth;
	}
	@Override
	public boolean hasNext() {
		ensureBuffer();
		return !buffer.isEmpty();
	}
	@Override
	public KmerPathNode next() {
		ensureBuffer();
		return buffer.poll();
	}
	private void ensureBuffer() {
		while (underlying.hasNext() && buffer.peek().startPosition() + maxSupportWidth >= underlying.peek().startPosition()) {
			process(underlying.next());
		}
		if (!underlying.hasNext()) {
			// flush everything
			flushEndingBefore(Integer.MAX_VALUE);
		}
	}
	/**
	 * flushes PathNodes that end before the given position
	 * @param position
	 */
	private void flushEndingBefore(int position) {
		throw new RuntimeException("NYI");
	}
	private void process(KmerAggregateNode next) {
		// complete all nodes that end before the start
		flushEndingBefore(next.startPosition());
		// Incorporate this node into the active paths
		
		// Invariant: start position of all path nodes < current node
		
		// add as successor
		
		// split existing node
		
		// add as predecessor
		
	}
	@Override
	public void remove() {
		throw new UnsupportedOperationException();
	}
}
