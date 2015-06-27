package au.edu.wehi.idsv.debruijn.positional;

import java.util.ArrayDeque;
import java.util.Iterator;
import java.util.PriorityQueue;

import com.google.common.collect.Iterators;
import com.google.common.collect.Ordering;
import com.google.common.collect.PeekingIterator;
import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;
import com.google.common.primitives.Ints;


/**
 * Calls optimal contigs from a positional de Bruijn graph
 * 
 * @author Daniel Cameron
 *
 */
public class BestNonReferenceContigCaller {
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
	private final PriorityQueue<Contig> called = new PriorityQueue<Contig>(1024, Contig.ByScoreDesc);
	private final int maxEvidenceWidth;
	private int inputPosition;
	public BestNonReferenceContigCaller(
			Iterator<KmerPathNode> it,
			int maxEvidenceWidth) {
		this.underlying = Iterators.peekingIterator(it);
		this.maxEvidenceWidth = maxEvidenceWidth;
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
		for (TraversalNode head = frontier.peekFrontier(); head != null && head.node.endPosition() < inputPosition; head = frontier.peekFrontier()) {
			visit(frontier.pollFrontier());
		}
	}
	private void queueForProcessing(KmerPathNode node) {
		if (!node.isReference()) {
			unprocessedStartNodes.add(node);
		}
	}
	/**
	 * Adds starting paths for the given node to the graph
	 * 
	 * A starting path is a position interval in which either
	 * a) no predecessors exist, or
	 * b) at least one reference kmer predecessor exists
	 * 
	 * @param node
	 */
	private void addStartingPaths(KmerPathNode node) {
		assert(!node.isReference());
		PeekingIterator<KmerPathNode> startIt = Iterators.peekingIterator(node.prev().iterator());
		int start = node.startPosition(0);
		final int scopeEnd = node.endPosition(0);
		int nonReferenceCount = 0;
		int referenceCount = 0;
		// TODO: is using a linear array faster?
		PriorityQueue<KmerPathNode> active = new PriorityQueue<KmerPathNode>(3, KmerNode.ByEndPosition);
		while (start <= scopeEnd) {
			// advance
			while (startIt.hasNext() && startIt.peek().startPosition() < start) {
				KmerPathNode n = startIt.next();
				if (n.isReference()) referenceCount++;
				else nonReferenceCount++;
				active.add(n);
			}
			while (!active.isEmpty() && active.peek().endPosition() + 1 < start) {
				KmerPathNode n = active.poll();
				if (n.isReference()) referenceCount--;
				else nonReferenceCount--;
			}
			int end = scopeEnd;
			if (startIt.hasNext()) {
				end = Math.min(end, startIt.peek().startPosition());
			}
			if (!active.isEmpty()) {
				end = Math.min(end, active.peek().endPosition() + 1);
			}
			if (referenceCount > 0) {
				// start of anchored path
				frontier.memoize(new TraversalNode(new KmerPathSubnode(node, start, end), ANCHORED_SCORE));
			} else if (referenceCount == 0 && nonReferenceCount == 0) {
				// start of unanchored path
				frontier.memoize(new TraversalNode(new KmerPathSubnode(node, start, end), 0));
			}
			start = end + 1;
		}
	}
	private void visit(TraversalNode ms) {
		assert(ms.node.endPosition() < inputPosition); // successors must be fully defined
		RangeSet<Integer> terminalAnchor = null;
		for (KmerPathSubnode sn : ms.node.next()) {
			if (!sn.node().isReference()) {
				frontier.memoize(new TraversalNode(ms, sn));
			} else {
				if (terminalAnchor == null) {
					terminalAnchor = TreeRangeSet.create();
				}
				terminalAnchor.add(Range.closed(ms.node.firstKmerStartPosition(), ms.node.firstKmerEndPosition()));
			}
		}
		if (terminalAnchor != null) {
			// path has reference successor = path is anchored to the reference here
			for (Range<Integer> rs : terminalAnchor.asRanges()) {
				called.add(new Contig(new TraversalNode(ms, rs.lowerEndpoint(), rs.upperEndpoint()), true));
			}	
		}
		for (Range<Integer> rs : ms.node.nextPathRangesOfDegree(0).asRanges()) {
			// path has no successors = end of path
			called.add(new Contig(new TraversalNode(ms, rs.lowerEndpoint(), rs.upperEndpoint()), false));
		}
	}
	public static class Contig {
		public Contig(TraversalNode node, boolean hasReferenceSuccessor) {
			this.node = node;
			this.score = node.score + (hasReferenceSuccessor ? ANCHORED_SCORE : 0);			
		}
		/**
		 * terminal contig node
		 */
		public final TraversalNode node;
		/**
		 * Final score for contig (including any anchor scoring bonus)
		 */
		public final int score;
		public ArrayDeque<KmerPathSubnode> toSubnodePath() {
			ArrayDeque<KmerPathSubnode> contigPath = new ArrayDeque<KmerPathSubnode>();
			KmerPathSubnode last = node.node;
			contigPath.add(last);
			for (TraversalNode n = node.prev; n != null; n = n.prev) {
				KmerPathSubnode current = n.node.givenNext(last);
				last = current;
				contigPath.addFirst(current);
			}
			return contigPath;
		}
		@Override
		public String toString() {
			return String.format("Path Score %d, %s", score, node);
		}
		public static final Ordering<Contig> ByScore = new Ordering<Contig>() {
			@Override
			public int compare(Contig left, Contig right) {
				return Ints.compare(left.score, right.score);
			}};
		public static final Ordering<Contig> ByScoreDesc = new Ordering<Contig>() {
			@Override
			public int compare(Contig left, Contig right) {
				return Ints.compare(right.score, left.score);
			}};
	}
	public ArrayDeque<KmerPathSubnode> bestContig() {
		while (underlying.hasNext() && (called.isEmpty() || frontier.peekFrontier() == null || called.peek().node.node.endPosition() > frontier.peekFrontier().node.firstKmerStartPosition() - maxEvidenceWidth)) {
			// if our best assembly could share evidence with an assembly that could be better, then we need to process more
			advance();
		}
		if (!underlying.hasNext()) {
			// final advance to end of input
			advance();
			assert(frontier.peekFrontier() == null);
		}
		Contig best = called.poll();
		if (best == null) return null;
		return best.toSubnodePath();
	}
}
