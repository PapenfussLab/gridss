package au.edu.wehi.idsv.debruijn.positional;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Stack;

import au.edu.wehi.idsv.AssemblyParameters;

/**
 * Simplifies the graph until the number of nodes used to represent the graph cannot be further reduced
 * 
 *  Without loss of generality, consider merging two nodes N_1, N_2 with N_1 \in next(M): A holds before merge
 *  - merging nodes N_1 and N_2 into N' will cause A to be violated iff:
 *   - adjacent node M can be collapsed into N', but could not be collapsed into N_1.
 * To merge nodes:
 *  start(N_1) = start(N_2)
 *  end(N_1) = end(N_2)
 *  kmers with hamming distance threshold
 *  |next(N_12)| >= |next(N_1)| degree does not decrease
 *  M degree reduced by 1 iff N_2 \in next(M)
 *  	N_12 can be merged with M 
 *   *  
 *  Consider the following graph:
 *  [11-15] 1 GTAC (node A)
 *  [16-20] 2 GTAC (node B)
 *  [11-15] 1 GGAC (node C)
 *  [10-19] 2 GGTA (node D)
 *  
 *  We merges C into A:
 *  [11-15] 2 GTAC (node A,C)
 *  [16-20] 2 GTAC (node B)
 *  [10-19] 2 GGTA (node D)
 *  
 *  Merging AC with B (since they differ only in validity intervals, and these intervals are adjacent) results in:  
 *  [11-20] 2 GTAC (node A,B,C)
 *  [10-19] 2 GGTA (node D)
 *  
 *  Merging ABC with D:
 *  [10-19] 4 GGTAC (nodes A,B,C,D)
 *  
 * A larger example would result in additional nodes being further merged.
 * The number of possible merges is unbounded, but since each merge results in a node containing N_1,
 * width of the merge is bounded by maxPathLength in either direction
 * 
 * @author cameron.d
 *
 */
public class PathSimplificationIterator implements Iterator<KmerPathNode> {
	private final HashMap<KmerNode.KmerStartPositionEqualityWrapper, KmerPathNode> startLookup = new HashMap<KmerNode.KmerStartPositionEqualityWrapper, KmerPathNode>();
	private final HashMap<KmerNode.KmerStartPositionEqualityWrapper, KmerPathNode> endLookup = new HashMap<KmerNode.KmerStartPositionEqualityWrapper, KmerPathNode>();
	private final int maxSupportWidth;
	private final int maxPathLength;
	private final Iterator<KmerPathNode> underlying;
	public PathSimplificationIterator(
			Iterator<KmerPathNode> it,
			AssemblyParameters ap,
			int maxSupportWidth) {
		this.underlying = it;
		this.maxSupportWidth = maxSupportWidth;
		this.maxPathLength = ap.positionalMaxPathLengthInBases(maxSupportWidth);
	}
	private void compressNeighbours(KmerPathNode root) {
		Stack<KmerPathNodeSnapshot> toProcess = new Stack<KmerPathNodeSnapshot>();
		toProcess.add(new KmerPathNodeSnapshot(root));
		while (!toProcess.isEmpty()) {
			KmerPathNodeSnapshot snapshot = toProcess.pop();
			if (snapshot.isUnchanged()) {
				KmerPathNode node = snapshot.node;
				if (node.next().size() == 1) {
					KmerPathNode candidate = node.next().get(0);
					if (candidate.prev().size() == 1
							&& candidate.startPosition(0) == node.endPosition() + 1
							&& candidate.endPosition(0) == node.endPosition() + 1
							&& node.length() + candidate.length() <= maxPathLength) {
						// we can concatenate these nodes together
						candidate.prepend(node);
						toProcess.add(new KmerPathNodeSnapshot(candidate));
						continue;
					}
				}
				if (node.prev().size() == 1) {
					KmerPathNode candidate = node.prev().get(0);
					if (candidate.next().size() == 1
							&& candidate.startPosition() + 1 == node.endPosition(0)
							&& candidate.endPosition() + 1 == node.endPosition(0)
							&& node.length() + candidate.length() <= maxPathLength) {
						// we can concatenate these nodes together
						node.prepend(candidate);
						toProcess.add(new KmerPathNodeSnapshot(node));
						continue;
					}
				}
			}
		}
	}

	@Override
	public boolean hasNext() {
		// TODO Auto-generated method stub
		throw new RuntimeException("NYI");
	}

	@Override
	public KmerPathNode next() {
		// TODO Auto-generated method stub
		throw new RuntimeException("NYI");
	}

	@Override
	public void remove() {
		throw new UnsupportedOperationException();
	}
}
