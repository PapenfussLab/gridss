package au.edu.wehi.idsv.debruijn.positional;

import au.edu.wehi.idsv.util.RangeUtil;
import com.google.common.collect.*;
import org.apache.commons.lang3.StringUtils;

import java.util.ArrayDeque;
import java.util.Iterator;
import java.util.List;

public abstract class KmerPathNodeBasePath {
	private final TraversalNode root;
	private final boolean traverseForward;
	private final int maxPathLength;
	protected TraversalNode rootNode() { return root; }
	public boolean traversingForward() { return traverseForward; }
	public KmerPathNodeBasePath(KmerPathSubnode node, int maxPathLength, boolean traverseForward) {
		this.traverseForward = traverseForward;
		this.maxPathLength = maxPathLength;
		this.root = new TraversalNode(null, node);
	}
	public class TraversalNode implements Iterable<TraversalNode> {
		public TraversalNode(TraversalNode parent, KmerPathSubnode node) {
			this.parent = parent;
			this.node = node;
			if (parent == null) {
				this.pathLength = node.length();
				this.pathWeight = node.weight();
			} else {
				this.pathLength = node.length() + parent.pathLength;
				this.pathWeight = node.weight() + parent.pathWeight;
			}
		}
		private final KmerPathSubnode node;
		private final TraversalNode parent;
		private final int pathLength;
		private final int pathWeight;
		public int pathLength() { return pathLength; }
		public int pathWeight() { return pathWeight; }
		public TraversalNode parent() { return parent; }
		public KmerPathSubnode node() { return node; }
		/**
		 * Returns the first kmer terminal ranges of this node 
		 * @return first kmer positions of this node that have no successors
		 */
		public RangeSet<Integer> terminalRanges() {
			return traversingForward() ? node.nextPathRangesOfDegree(KmerPathSubnode.NO_EDGES) : node.prevPathRangesOfDegree(KmerPathSubnode.NO_EDGES);
		}
		/**
		 * Finds all the terminal leaves that are subsets of this path
		 * 
		 * For each node, the set of leaf intervals is the set of intervals
		 * for which the only predecessor or successors are nodes on 
		 * our traversal path
		 * Consider the path ABC
		 * A -> B <- C           this is a leaf interval
		 * A -> B <- C -> D      CD makes this not a leaf interval
		 * A -> B <- C           BE makes this not a leaf interval
		 * A -> B \
		 * A       E              
		 * A
		 *  
		 * @return intervals for which this path is a terminal leaf. Intervals are for position of the path anchor kmer
		 */
		public RangeSet<Integer> terminalLeafAnchorRanges() {
			RangeSet<Integer> ranges = toAnchorPosition(terminalRanges());
			TraversalNode tn = this;
			KmerPathSubnode sn = node;
			while (tn != null && !ranges.isEmpty()) {
				if (tn != this) {
					// only a single successor
					RangeSet<Integer> singleSuccessor = traversingForward() ? sn.nextPathRangesOfDegree(KmerPathSubnode.SINGLE_EDGE) : sn.prevPathRangesOfDegree(KmerPathSubnode.SINGLE_EDGE);
					singleSuccessor = tn.toAnchorPosition(singleSuccessor);
					ranges = RangeUtil.intersect(ranges, singleSuccessor);
				}
				// only a single ancestor
				RangeSet<Integer> singleAncestor = traversingForward() ? sn.prevPathRangesOfDegree(KmerPathSubnode.SINGLE_EDGE) : sn.nextPathRangesOfDegree(KmerPathSubnode.SINGLE_EDGE);
				singleAncestor = tn.toAnchorPosition(singleAncestor);
				ranges = RangeUtil.intersect(ranges, singleAncestor);
				tn = tn.parent;
				if (tn != null) {
					sn = traversingForward() ? tn.node.givenNext(sn) :  tn.node.givenPrev(sn);
				}
			}
			return ranges;
		}
		public int startPositionOfAnchorKmer() {
			return firstNodeKmerToAnchorPosition(node().firstStart()); 
		}
		public int endPositionOfAnchorKmer() {
			return firstNodeKmerToAnchorPosition(node().firstEnd());
		}
		/**
		 * Translates a path position of the first kmer of this to the corresponding anchor kmer
		 * position  
		 * @param pos node first kmer position
		 * @return position of the first kmer encountered in path traversal 
		 */
		private int firstNodeKmerToAnchorPosition(int pos) {
			if (traversingForward()) {
				return pos + node.length() - pathLength;
			} else {
				return pos + pathLength - 1;
			}
		}
		private int anchorTofirstNodeKmerPosition(int pos) {
			return pos - firstNodeKmerToAnchorPosition(0);
		}
		/**
		 * Translates a path position of the first kmer of this to the corresponding
		 * position of the first kmer of the path root  
		 * @param pos node first kmer position
		 * @return root first kmer position 
		 */
		private RangeSet<Integer> toAnchorPosition(RangeSet<Integer> rs) {
			RangeSet<Integer> rootRanges = TreeRangeSet.create();
			for (Range<Integer> r : rs.asRanges()) {
				rootRanges.add(Range.closed(firstNodeKmerToAnchorPosition(r.lowerEndpoint()), firstNodeKmerToAnchorPosition(r.upperEndpoint())));
			}
			return rootRanges;
		}
		public ArrayDeque<KmerPathSubnode> asSubnodes() {
			ArrayDeque<KmerPathSubnode> queue = new ArrayDeque<KmerPathSubnode>();
			TraversalNode tn = this;
			KmerPathSubnode sn = node;
			while (tn != null) {
				if (traversingForward()) {
					queue.addFirst(sn);
				} else {
					queue.addLast(sn);
				}
				tn = tn.parent;
				if (tn != null) {
					sn = traversingForward() ? tn.node.givenNext(sn) :  tn.node.givenPrev(sn);
				}
			}
			return queue;
		}
		public TraversalNode overlapping(TraversalNode path) {
			int startAnchor = Math.max(startPositionOfAnchorKmer(), path.startPositionOfAnchorKmer());
			int endAnchor = Math.min(endPositionOfAnchorKmer(), path.endPositionOfAnchorKmer());
			if (startAnchor > endAnchor) return null; // no overlap
			if (startAnchor == startPositionOfAnchorKmer() && endAnchor == endPositionOfAnchorKmer()) {
				return this;
			} else {
				return new TraversalNode(parent, new KmerPathSubnode(node().node(), anchorTofirstNodeKmerPosition(startAnchor), anchorTofirstNodeKmerPosition(endAnchor)));
			}
		}
		public TraversalNode firstTerminalLeaf() {
			RangeSet<Integer> rs = terminalLeafAnchorRanges();
			if (rs == null || rs.isEmpty()) return null;
			Range<Integer> firstRange = rs.asRanges().iterator().next();
			int start = firstRange.lowerEndpoint();
			int end = firstRange.upperEndpoint();
			start = anchorTofirstNodeKmerPosition(start);
			end = anchorTofirstNodeKmerPosition(end);
			return new TraversalNode(parent, new KmerPathSubnode(node().node(), start, end));
		}
		@Override
		public Iterator<TraversalNode> iterator() {
			List<KmerPathSubnode> adj = traversingForward() ? node.next() : node.prev();
			return new TraversalNodeIterator(this, adj.iterator());
		}
		public String toString() {
			return StringUtils.stripEnd(asSubnodes().toString().replace(",", "\n").substring(1), "]");
		}
	}
	private class TraversalNodeIterator implements Iterator<TraversalNode> {
		private final TraversalNode node;
		private PeekingIterator<KmerPathSubnode> successors;
		public TraversalNodeIterator(TraversalNode node, Iterator<KmerPathSubnode> successors) {
			this.node = node;
			this.successors = Iterators.peekingIterator(successors);
		}
		private void ensureNext() {
			while (successors.hasNext() && successors.peek().node().length() + node.pathLength > maxPathLength) {
				// skip paths longer than our maximum path length
				successors.next();
			}
		}
		@Override
		public boolean hasNext() {
			ensureNext();
			return successors.hasNext();
		}
		@Override
		public TraversalNode next() {
			ensureNext();
			return new TraversalNode(node, successors.next());
		}
		@Override
		public void remove() {
			throw new UnsupportedOperationException();
		}
	}
}
