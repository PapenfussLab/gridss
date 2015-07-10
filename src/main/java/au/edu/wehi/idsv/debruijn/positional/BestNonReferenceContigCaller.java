package au.edu.wehi.idsv.debruijn.positional;

import java.util.ArrayDeque;
import java.util.Iterator;
import java.util.PriorityQueue;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Iterators;
import com.google.common.collect.Ordering;
import com.google.common.collect.PeekingIterator;
import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;


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
	private final PriorityQueue<KmerPathNode> unprocessedStartNodes = new PriorityQueue<KmerPathNode>(KmerNodeUtil.ByLastEnd);
	/**
	 * Assembled contigs. The best contig is called once all potential alternate contigs
	 * involving the evidence used to construct the contig have also been assembled.
	 */
	private final PriorityQueue<Contig> called = new PriorityQueue<Contig>(Contig.ByScoreDescPosition);
	private final int maxEvidenceWidth;
	private int inputPosition;
	private long consumed = 0;
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
			inputPosition = underlying.peek().firstStart();
		}
		advanceUnderlying();
		advanceUnprocessed();
		advanceFrontier();
	}
	/**
	 * Loads records from the underlying stream up to and including the current inputPosition.
	 */
	private void advanceUnderlying() {
		while (underlying.hasNext() && underlying.peek().firstStart() <= inputPosition) {			
			KmerPathNode nextRecord = underlying.next();
			consumed++;
			queueForProcessing(nextRecord);
		}
	}
	private void advanceUnprocessed() {
		// final kmer ending before inputPosition means that all adjacent nodes have been loaded
		while (!unprocessedStartNodes.isEmpty() && unprocessedStartNodes.peek().lastEnd() < inputPosition) {
			KmerPathNode unprocessed = unprocessedStartNodes.poll();
			addStartingPaths(unprocessed);
		}
	}
	private void advanceFrontier() {
		for (TraversalNode head = frontier.peekFrontier(); head != null && head.node.lastEnd() < inputPosition; head = frontier.peekFrontier()) {
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
		int start = node.firstStart();
		final int scopeEnd = node.firstEnd();
		int nonReferenceCount = 0;
		int referenceCount = 0;
		// TODO: is using a linear array faster?
		PriorityQueue<KmerPathNode> active = new PriorityQueue<KmerPathNode>(3, KmerNodeUtil.ByLastEnd);
		while (start <= scopeEnd) {
			// advance
			while (startIt.hasNext() && startIt.peek().lastStart() < start) {
				KmerPathNode n = startIt.next();
				if (n.isReference()) referenceCount++;
				else nonReferenceCount++;
				active.add(n);
			}
			while (!active.isEmpty() && active.peek().lastEnd() + 1 < start) {
				KmerPathNode n = active.poll();
				if (n.isReference()) referenceCount--;
				else nonReferenceCount--;
			}
			int end = scopeEnd;
			if (startIt.hasNext()) {
				end = Math.min(end, startIt.peek().lastStart());
			}
			if (!active.isEmpty()) {
				end = Math.min(end, active.peek().lastEnd() + 1);
			}
			if (referenceCount > 0) {
				// start of anchored path
				frontier.memoize(new MemoizedTraversalNode(new KmerPathSubnode(node, start, end), ANCHORED_SCORE));
			} else if (referenceCount == 0 && nonReferenceCount == 0) {
				// start of unanchored path
				frontier.memoize(new MemoizedTraversalNode(new KmerPathSubnode(node, start, end), 0));
			}
			start = end + 1;
		}
	}
	private void visit(TraversalNode ms) {
		assert(ms.node.lastEnd() < inputPosition); // successors must be fully defined
		RangeSet<Integer> terminalAnchor = null;
		for (KmerPathSubnode sn : ms.node.next()) {
			if (!sn.node().isReference()) {
				frontier.memoize(new MemoizedTraversalNode(ms, sn));
			} else {
				if (terminalAnchor == null) {
					terminalAnchor = TreeRangeSet.create();
				}
				terminalAnchor.add(Range.closed(ms.node.firstStart(), ms.node.firstEnd()));
			}
		}
		if (terminalAnchor != null) {
			// path has reference successor = path is anchored to the reference here
			for (Range<Integer> rs : terminalAnchor.asRanges()) {
				called.add(new Contig(new TraversalNode(ms, rs.lowerEndpoint(), rs.upperEndpoint()), true));
			}	
		}
		for (Range<Integer> rs : ms.node.nextPathRangesOfDegree(KmerPathSubnode.NO_EDGES).asRanges()) {
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
			return node.toSubnodeNextPath();
		}
		@Override
		public String toString() {
			return String.format("Path Score %d, %s", score, node);
		}
		public static final Ordering<Contig> ByScoreDescPosition = new Ordering<Contig>() {
			@Override
			public int compare(Contig left, Contig right) {
				return ComparisonChain.start()
						.compare(right.score, left.score)
						.compare(left.node.node.lastStart(), right.node.node.lastStart())
						.result();
			}};
	}
	/**
	 * Determines whether the contig could overlap evidence with another contig
	 * that we have not yet generated. We need to ensure that the contig is
	 * at least maxEvidenceWidth away from:
	 * a) all partially generated contig
	 * b) all unprocessed contig start locations
	 * c) any additional nodes not ryet read from the underlying stream
	 *  
	 * @param contig
	 * @return true if the contig will not share evidence with future contigs, false otherwise
	 */
	private boolean contigDoesNotShareEvidenceWithUnprocessed(Contig contig) {
		assert(contig != null);
		int contigLastEnd = contig.node.node.lastEnd();
		int frontierFirstStart = frontier.peekFrontier() == null ? Integer.MAX_VALUE : frontier.peekFrontier().node.firstStart();
		int unprocessedFirstNodeLastEnd = unprocessedStartNodes.isEmpty() ? Integer.MAX_VALUE : unprocessedStartNodes.peek().lastEnd();
		int unprocessedFirstStart = unprocessedFirstNodeLastEnd - maxEvidenceWidth; // node could contain entire evidence
		int firstStart = Math.min(frontierFirstStart, Math.min(unprocessedFirstStart, inputPosition));
		return contigLastEnd < firstStart - maxEvidenceWidth; // evidence could overlap just contig last end
	}
	public ArrayDeque<KmerPathSubnode> bestContig() {
		// FIXME: add hard safety bounds to the width of the loaded graph
		// since the size is technically unbounded
		while (underlying.hasNext() && (called.isEmpty() || !contigDoesNotShareEvidenceWithUnprocessed(called.peek()))) {
			advance();
		}
		if (!underlying.hasNext()) {
			// final advance to end of input
			advance();
			assert(frontier.peekFrontier() == null);
		}
		Contig best = called.poll();
		if (best == null) {
			assert(!underlying.hasNext());
			return null;
		}
		return best.toSubnodePath();
	}
	public int tracking_contigCount() {
		return called.size();
	}
	public int tracking_contigFirstPosition() {
		return called.stream().mapToInt(c -> c.toSubnodePath().getFirst().firstStart()).min().orElse(Integer.MAX_VALUE);
	}
	public long tracking_underlyingConsumed() {
		return consumed;
	}
	public int tracking_memoizedNodeCount() {
		return frontier.tracking_memoizedNodeCount();
	}
	public int tracking_frontierSize() {
		return frontier.tracking_frontierSize();
	}
	public int tracking_unprocessedStartNodeCount() {
		return unprocessedStartNodes.size();
	}
	
}
