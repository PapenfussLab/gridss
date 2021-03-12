package au.edu.wehi.idsv.debruijn.positional;

import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures.KmerNodeByFirstStartPriorityQueue;
import au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures.KmerNodeByLastEndPriorityQueue;
import au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures.KmerNodeNonOverlappingLookup;
import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;
import htsjdk.samtools.util.Log;

import java.util.*;

/**
 * Transforms a start position sorted sequence of non-overlapping KmerAggregateNode to a
 * start position sorted sequence of KmerPathNode with graph edges 
 * @author Daniel Cameron
 *
 */
public class PathNodeIterator implements Iterator<KmerPathNode> {
	private static final Log log = Log.getInstance(PathNodeIterator.class);
	private final PeekingIterator<? extends KmerNode> underlying;
	private final int maxNodeLength;
	private final int k;
	private final KmerNodeNonOverlappingLookup<KmerNode> edgeLookup;
	/**
	 * Nodes that have not yet had all edges defined.
	 * 
	 * Edges are fully defined when the inputPosition is greater than the end position of the node
	 * - latest first kmer start position of node is equal to its endPosition
	 * - prev nodes must end before current not first kmer end position
	 * - next nodes must start no later than our end position + 1
	 */
	private final Queue<KmerNode> activeNodes = Defaults.USE_OPTIMISED_ASSEMBLY_DATA_STRUCTURES ? new KmerNodeByLastEndPriorityQueue<>(16) : new PriorityQueue<>(KmerNodeUtil.ByLastEnd);
	private final Queue<KmerPathNode> pathNodes = Defaults.USE_OPTIMISED_ASSEMBLY_DATA_STRUCTURES ? new KmerNodeByFirstStartPriorityQueue<>(16) : new PriorityQueue<>(1024, KmerNodeUtil.ByFirstStart);
	private int inputPosition = Integer.MIN_VALUE;
	/**
	 * Maximum width of a single node. This is calculated from the input sequence
	 */
	private int maxNodeWidth = 0;
	private long consumed = 0;
	public PathNodeIterator(Iterator<? extends KmerNode> it, int maxPathLength, int k) {
		if (maxPathLength < 1) throw new IllegalArgumentException("Path length must be positive");
		this.underlying = Iterators.peekingIterator(it);
		this.maxNodeLength = maxPathLength;
		this.k = k;
		this.edgeLookup = new KmerNodeNonOverlappingLookup<>(k);
	}
	@Override
	public boolean hasNext() {
		boolean hasMore = !pathNodes.isEmpty() || !activeNodes.isEmpty() || underlying.hasNext();
		if (!hasMore) {
			assert(pathNodes.isEmpty());
			assert(activeNodes.isEmpty());
			assert(!underlying.hasNext());
			assert(edgeLookup.isEmpty());
		}
		return hasMore;
	}
	private void advance() {
		inputPosition = underlying.peek().lastStart();
		while (underlying.hasNext() && underlying.peek().lastStart() == inputPosition) {
			KmerNode node = underlying.next();
			consumed++;
			edgeLookup.add(node);
			activeNodes.add(node);
			maxNodeWidth = Math.max(maxNodeWidth, node.width());
		}
		if (Defaults.SANITY_CHECK_ASSEMBLY_GRAPH) {
			sanityCheck();
		}
	}
	private void merge() {
		while (!activeNodes.isEmpty() && activeNodes.peek().lastEnd() < inputPosition) {
			KmerNode node = activeNodes.poll();
			merge(node);
		}
		if (Defaults.SANITY_CHECK_ASSEMBLY_GRAPH) {
			sanityCheck();
		}
	}
	private void merge(KmerNode right) {
		// can we merge into an earlier KmerNode?
		KmerNode left = edgeLookup.getUniqueFullWidthPredecessor(right);
		if (left != null) {
			//assert(edgeLookup.prevNodes(right).size() == 1); // TEMP HACK
			assert (left instanceof KmerPathNode); // must have already processed our previous node
			if (left.length() < maxNodeLength &&
					right.isReference() == left.isReference() &&
					edgeLookup.getUniqueFullWidthSuccessor(left) == right) {
				//assert(edgeLookup.nextNodes(left).size() == 1); // TEMP HACK
				assert(KmerEncodingHelper.isNext(k, left.lastKmer(), right.firstKmer()));
				KmerPathNode pn = (KmerPathNode) left;
				edgeLookup.adjustForMerge(left, right);
				pn.append(right);
				return;
			}
		}
		// couldn't merge into a previous path = new path
		KmerPathNode pn = new KmerPathNode(right);
		pathNodes.add(pn);
		edgeLookup.replace(right, pn);
	}
	@Override
	public KmerPathNode next() {
		while (underlying.hasNext() && (pathNodes.isEmpty() || !edgesCanBeFullyDefined(pathNodes.peek()))) {
			advance();
			merge();
		}
		if (!underlying.hasNext()) {
			// finish off
			inputPosition = Integer.MAX_VALUE;
			merge();
		}
		KmerPathNode pn = pathNodes.poll();
		assert(pn != null);
		assert(edgesCanBeFullyDefined(pn));
		// populate remaining edges
		for (KmerNode n : edgeLookup.prevNodes(pn)) {
			assert(n instanceof KmerPathNode);
			KmerPathNode.addEdge((KmerPathNode)n, pn);
		}
		for (KmerNode n : edgeLookup.nextNodes(pn)) {
			assert(n instanceof KmerPathNode);
			if (n != pn) {
				KmerPathNode.addEdge(pn, (KmerPathNode)n);
			}
		}
		edgeLookup.remove(pn);
		assert(pn.length() <= maxNodeLength);
		if (Defaults.SANITY_CHECK_ASSEMBLY_GRAPH) {
			sanityCheck();
		}
		return pn;
	}
	private boolean edgesCanBeFullyDefined(KmerPathNode node) {
		// to still be unprocessed, the node has to end at or after inputPosition
		// therefore, earliest start position of a unprocessed node
		// is based on the maximum width of a node
		int earliestFirstKmerStartPositionOfActiveNode = inputPosition - maxNodeWidth;
		return node.lastEnd() < earliestFirstKmerStartPositionOfActiveNode;
	}
	@Override
	public void remove() {
		throw new UnsupportedOperationException();
	}
	public boolean sanityCheck() {
		edgeLookup.sanityCheck();
		for (KmerNode n : activeNodes) {
			assert(!pathNodes.contains(n));
		}
		for (KmerPathNode n : pathNodes) {
			assert(!activeNodes.contains(n));
			n.sanityCheck();
		}
		assert(activeNodes.stream().distinct().count() == activeNodes.size());
		assert(pathNodes.stream().distinct().count() == pathNodes.size());
		return true;
	}
	public int tracking_processedSize() {
		return pathNodes.size();
	}
	public int tracking_activeSize() {
		return activeNodes.size();
	}
	public int tracking_inputPosition() {
		return inputPosition;
	}
	public long tracking_underlyingConsumed() {
		return consumed;
	}
	public int tracking_edgeLookupSize() {
		return edgeLookup.size();
	}
	public int tracking_pathNodeEdgeLookupSize() {
		return -1;
	}
	public int tracking_edgeLookupMaxKmerNodeCount() {
		return -1;
	}
	public int tracking_pathNodeEdgeLookupMaxKmerNodeCount() {
		return -1;
	}
}
