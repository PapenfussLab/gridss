package au.edu.wehi.idsv.debruijn.positional;

import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap;

import java.util.Iterator;
import java.util.PriorityQueue;

import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.debruijn.positional.AggregateNodeIterator.KmerNodeAggregator.KmerNodeAggregatorSnapshot;

import com.google.common.collect.Iterators;
import com.google.common.collect.Ordering;
import com.google.common.collect.PeekingIterator;
import com.google.common.primitives.Ints;
import com.google.common.primitives.Longs;

/**
 * Transforms a start position sorted sequence of KmerNodes to a
 * start position sorted sequence of non-overlapping KmerAggregateNodes
 * @author Daniel Cameron
 *
 */
public class AggregateNodeIterator implements Iterator<KmerNode> {
	private final PeekingIterator<? extends KmerNode> underlying;
	private PriorityQueue<ImmutableKmerNode> outputSortBuffer = new PriorityQueue<ImmutableKmerNode>(1024, KmerNodeUtil.ByStartPosition);
	private Long2ObjectOpenHashMap<KmerNodeAggregator> byKmer = new Long2ObjectOpenHashMap<KmerNodeAggregator>();
	private PriorityQueue<KmerNodeAggregatorSnapshot> byStartPosition = new PriorityQueue<KmerNodeAggregatorSnapshot>(1024, BySnapshotEndPosition);
	public AggregateNodeIterator(Iterator<? extends KmerNode> it) {
		this.underlying = Iterators.peekingIterator(it);
	}
	@Override
	public boolean hasNext() {
		ensureBuffer();
		return !outputSortBuffer.isEmpty();
	}
	@Override
	public KmerNode next() {
		ensureBuffer();
		return outputSortBuffer.poll();
	}
	private void ensureBuffer() {
		// we can emit whenever there are no unprocessed or incomplete intervals
		// before our current interval
		while (underlying.hasNext() && (outputSortBuffer.isEmpty() ||
				outputSortBuffer.peek().startPosition() >= underlying.peek().startPosition() ||
				(!byStartPosition.isEmpty() && outputSortBuffer.peek().endPosition() >= byStartPosition.peek().snapshotEnd))) {
			process(underlying.peek().startPosition());
		}
		if (!underlying.hasNext()) {
			// flush everything
			flushBefore(Integer.MAX_VALUE);
		} else {
			flushBefore(underlying.peek().startPosition());
		}
		if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
			assert(sanityCheck());
		}
	}
	private void process(int position) {
		while (underlying.hasNext() && underlying.peek().startPosition() <= position) {
			KmerNode n = underlying.next();
			assert(position == Integer.MAX_VALUE || n.startPosition() == position); // input should be sorted by start position
			long kmer = n.kmer();
			KmerNodeAggregator ag = byKmer.get(kmer);
			if (ag == null) {
				ag = new KmerNodeAggregator(kmer);
				byKmer.put(kmer, ag);
			}
			ag.add(n);
			byStartPosition.add(ag.new KmerNodeAggregatorSnapshot());
		}
	}
	/**
	 * Flush all aggregate nodes starting before the given position
	 * @param position
	 */
	private void flushBefore(int position) {
		byStartPosition_removeInvalidHead();
		while (!byStartPosition.isEmpty() && byStartPosition.peek().snapshotEnd < position) {
			KmerNodeAggregatorSnapshot snapshot = byStartPosition.poll();
			assert(snapshot.isValid());
			KmerNodeAggregator ag = snapshot.aggregator();
			ag.advanceTo(position - 1);
			if (ag.isEmpty()) {
				byKmer.remove(ag.kmer);
			} else {
				KmerNodeAggregatorSnapshot newSnapshot = ag.new KmerNodeAggregatorSnapshot();
				assert(newSnapshot.snapshotEnd >= position);
				byStartPosition.add(newSnapshot);
			}
			byStartPosition_removeInvalidHead();
		}
	}
	private static final Ordering<KmerNodeAggregatorSnapshot> BySnapshotEndPosition = new Ordering<KmerNodeAggregatorSnapshot>() {
		@Override
		public int compare(KmerNodeAggregatorSnapshot left, KmerNodeAggregatorSnapshot right) {
			return Ints.compare(left.snapshotEnd, right.snapshotEnd);
		}
	};
	/**
	 * Generates KmerAggregateNode from an underlying sequence of KmerSupportNodes in ascending starting position
	 * 
	 * @author Daniel Cameron
	 *
	 */
	class KmerNodeAggregator implements Comparable<KmerNodeAggregator> {
		public class KmerNodeAggregatorSnapshot {
			public KmerNodeAggregatorSnapshot() {
				this.snapshotEnd = KmerNodeAggregator.this.end();
			}
			public final int snapshotEnd;
			/**
			 * Determines whether the snapshot is still valid
			 * @return
			 */
			public boolean isValid() {
				return this.snapshotEnd == KmerNodeAggregator.this.end();
			}
			// what's the correct syntax for doing this from outside the class?
			public KmerNodeAggregator aggregator() { return KmerNodeAggregator.this; }
		}
		public KmerNodeAggregator(long kmer) {
			this.kmer = kmer;
		}
		public boolean isEmpty() { return active.isEmpty(); }
		/**
		 * KmerNodes in the currently active aggregation interval
		 */
		private PriorityQueue<KmerNode> active = new PriorityQueue<KmerNode>(8, KmerNodeUtil.ByEndPosition);
		/**
		 * Start position of currently active aggregation interval
		 */
		private int start = Integer.MIN_VALUE;
		/**
		 * Weight of currently active aggregation interval
		 */
		private int weight = 0;
		/**
		 * Number of active reference KmerNode 
		 */
		private int referenceCount = 0;
		/**
		 * Advances to the next node, adding aggregate nodes to the given collection
		 * @param node next node
		 * @param emitTo collection to emit aggregate records to
		 */
		private final long kmer;
		public int end() {
			if (active.isEmpty()) return Integer.MAX_VALUE;
			return active.peek().endPosition();
		}
		public void add(KmerNode node) {
			assert(node.kmer() == kmer);
			assert(node.startPosition() >= start);
			advanceTo(node.startPosition() - 1);
			if (weight > 0 && start < node.startPosition()) {
				outputSortBuffer.add(new ImmutableKmerNode(kmer, start, node.startPosition() - 1, referenceCount > 0, weight));
			}
			start = node.startPosition();
			if (node.isReference()) {
				referenceCount++;
			}
			weight += node.weight();
			active.add(node);
		}
		/**
		 * Process up to and including the given position
		 * @param position final processing position
		 * @param emitTo collection to emit aggregate records to
		 */
		public void advanceTo(int position) {
			while (!active.isEmpty() && active.peek().endPosition() <= position) {
				int end = active.peek().endPosition();
				outputSortBuffer.add(new ImmutableKmerNode(kmer, start, end, referenceCount > 0, weight));
				while (!active.isEmpty() && active.peek().endPosition() == end) {
					KmerNode endingHere = active.poll();
					assert(endingHere.kmer() == kmer);
					weight -= endingHere.weight();
					if (endingHere.isReference()) {
						referenceCount--;
					}
				}
				start = end + 1;
			}
		}
		@Override
		public int compareTo(KmerNodeAggregator right) {
			return Longs.compare(kmer, right.kmer);
		}
	}
	@Override
	public void remove() {
		throw new UnsupportedOperationException();
	}
	private void byStartPosition_removeInvalidHead() {
		while (!byStartPosition.isEmpty() && !byStartPosition.peek().isValid()) {
			byStartPosition.poll();
		}
	}
	private boolean sanityCheck() {
		// kmer lookup is correct
		assert(byKmer.entrySet().stream().allMatch(kvp -> kvp.getKey() == kvp.getValue().kmer));
		// empty aggregators have been removed
		assert(byKmer.values().stream().allMatch(ag -> !ag.active.isEmpty()));
		// could have many start position entries, but only one position is valid (and even that could have duplicate entries)
		assert(byStartPosition.size() >= byKmer.size());
		assert(byStartPosition.stream().allMatch(snapshot -> !snapshot.isValid() || byKmer.containsKey(snapshot.aggregator().kmer)));
		if (outputSortBuffer.isEmpty()) {
			assert(byKmer.isEmpty());
			assert(byStartPosition.isEmpty());
			assert(!underlying.hasNext());
		}
		if (!byStartPosition.isEmpty()) {
			assert(byStartPosition.peek().isValid());
			if (underlying.hasNext()) {
				assert(byStartPosition.peek().snapshotEnd >= underlying.peek().startPosition());
			}
		}
		
		return true;
	}
}
