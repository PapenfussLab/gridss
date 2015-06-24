package au.edu.wehi.idsv.debruijn.positional;

import java.util.Iterator;
import java.util.PriorityQueue;

import au.edu.wehi.idsv.debruijn.positional.MemoizedTraverse.MemoizedNode;

import com.google.common.collect.Iterators;
import com.google.common.collect.Ordering;
import com.google.common.collect.PeekingIterator;
import com.google.common.collect.Range;
import com.google.common.primitives.Ints;


/**
 * Calls optimal contigs from a positional de Bruijn graph
 * 
 * @author Daniel Cameron
 *
 */
public class NonReferenceContigCaller {
	/**
	 * Since reference kmers are not scored, calculating 
	 * highest weighted results in a preference for paths
	 * ending at a RP with sequencing errors over a path
	 * anchored to the reference. 
	 * 
	 * To ensure that the anchored paths is score higher
	 * than the unanchored paths, paths anchored to the
	 * reference are given a score adjustment larger than
	 * the largest expected score.
	 */
	private static final int ANCHORED_SCORE = Integer.MAX_VALUE >> 2; 
	private final PeekingIterator<KmerPathNode> underlying;
	private final MemoizedTraverse frontier = new MemoizedTraverse();
	/**
	 * Potential starting nodes.
	 * We need to wait until all previous nodes are defined before checking
	 * for starting node intervals.
	 */
	private final PriorityQueue<KmerPathNode> unprocessedStartNodes = new PriorityQueue<KmerPathNode>(1024, KmerPathNode.ByEndPosition);
	private final PriorityQueue<Contig> called = new PriorityQueue<Contig>(1024, Contig.ByScore.reverse());
	private final int maxEvidenceWidth;
	private final int maxAnchorLength;
	private int inputPosition;
	public NonReferenceContigCaller(
			Iterator<KmerPathNode> it,
			int maxEvidenceWidth,
			int maxAnchorLength) {
		this.underlying = Iterators.peekingIterator(it);
		this.maxEvidenceWidth = maxEvidenceWidth;
		this.maxAnchorLength = maxAnchorLength;
	}
	public Contig next() {
		while (underlying.hasNext() && (called.isEmpty() || frontier.peekFrontier() == null || called.peek().node.node.endPosition() > frontier.peekFrontier().node.firstKmerStartPosition() - maxEvidenceWidth)) {
			// if our best assembly could share evidence with an assembly that could be better, then we need to process more
			advance();
		}
		Contig best = called.poll();
		// TODO: extend into anchor sequence in both directions
		// either greedy or a full MemoizedTraverse
		return best;
	}
	private void advance() {
		if (!underlying.hasNext()) {
			inputPosition = Integer.MAX_VALUE;
		} else {
			inputPosition = underlying.peek().startPosition(0);
		}
		advanceUnderlying();
		advanceUnprocessed();
		advanceFrontier();
	}
	/**
	 * Loads records from the underlying stream up to and including the current inputPosition.
	 */
	private void advanceUnderlying() {
		while (underlying.hasNext() && underlying.peek().startPosition(0) <= inputPosition) {			
			KmerPathNode nextRecord = underlying.next();
			queueForProcessing(nextRecord);
		}
	}
	private void advanceUnprocessed() {
		// final kmer ending before inputPosition means that all adjacent nodes have been loaded
		while (!unprocessedStartNodes.isEmpty() && unprocessedStartNodes.peek().endPosition() < inputPosition) {
			KmerPathNode unprocessed = unprocessedStartNodes.poll();
			addStartingPaths(unprocessed);
		}
	}
	private void advanceFrontier() {
		MemoizedNode head = frontier.peekFrontier();
		while (head != null && head.node.endPosition() < inputPosition) {
			visit(frontier.pollFrontier());
		}
	}
	private void queueForProcessing(KmerPathNode node) {
		if (!node.isReference()) {
			unprocessedStartNodes.add(node);
		}
	}
	private void addStartingPaths(KmerPathNode node) {
		assert(!node.isReference());
		KmerPathSubnode sn = new KmerPathSubnode(node);
		// find nodes that that have
		// a) only reference predecessors
		// b) no predecessors
		// add to memoization cache
		// add to frontier
	}
	private void visit(MemoizedNode ms) {
		assert(ms.node.endPosition() < inputPosition); // successors must be fully defined
		for (KmerPathSubnode sn : ms.node.next()) {
			if (!sn.node().isReference()) {
				frontier.memoize(new MemoizedNode(ms, sn));
			} else {
				// terminal node
				// TODO: handle overlapping/multiple
				called.add(new Contig(new MemoizedNode(ms, sn.firstKmerStartPosition(), sn.firstKmerEndPosition()), true));
			}
		}
		for (Range<Integer> rs : ms.node.nextPathRangesOfDegree(0).asRanges()) {
			called.add(new Contig(new MemoizedNode(ms, rs.lowerEndpoint(), rs.upperEndpoint()), false));
		}
	}
	public static class Contig {
		public Contig(MemoizedNode node, boolean hasReferenceSuccessor) {
			this.node = node;
			this.score = node.score + (hasReferenceSuccessor ? ANCHORED_SCORE : 0);			
		}
		/**
		 * terminal contig node
		 */
		public final MemoizedNode node;
		/**
		 * Final score for contig (including any anchor scoring bonus)
		 */
		public final int score;
		public static final Ordering<Contig> ByScore = new Ordering<Contig>() {
			@Override
			public int compare(Contig left, Contig right) {
				return Ints.compare(left.score, right.score);
			}};
	}
}
