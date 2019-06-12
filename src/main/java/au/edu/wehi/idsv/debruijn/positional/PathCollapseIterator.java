package au.edu.wehi.idsv.debruijn.positional;

import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.debruijn.DeBruijnGraph;
import au.edu.wehi.idsv.debruijn.DeBruijnSequenceGraphNodeUtil;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.debruijn.positional.KmerPathNodeBasePath.TraversalNode;
import au.edu.wehi.idsv.util.IntervalUtil;
import au.edu.wehi.idsv.util.SequenceUtil;

import java.util.*;

/**
 * Graph simplifier that merging similar paths
 * 
 * WARNING: this performs full tree to tree comparison of all children and
 * runs in exponential time
 * 
 * Input: KmerPathNodes in ascending order of start position of the first kmer
 * 
 * Output: KmerPathNodes in ascending order of start position of the first kmer
 * after graph simplification. Each output node is guaranteed to not be modified
 * by further graph reduction although adjacent node may be modified causing a
 * change in the node edge list after emission.
 *  
 * All input nodes are guaranteed to have all edges defined,
 * but this is not transitive: nodes that have not been returned by the input
 * are not guaranteed to either a) have all edges defined, or b) be of full length
 *  
 * When collapsing 2 path together, all nodes along both paths must be fully defined
 * so the edge list of adjacent nodes can be updated if the node is split.
 * 
 * Let G(pos) be the graph of all nodes such that first_kmer_start_position(n) <= pos
 * 
 * For branch collapse:
 * ... -* - * <- root node we traverse backwards looking for similar paths
 *        /
 * ...   *
 * 
 * Let n be a node such that next(n) = root
 * To collapse a path onto n, n needs to be split along:
 * a) length: alternate path may have nodes of shorter length,
 *     this requires this node to be split into shorter nodes
 * b) start/end: the alternate path may have a smaller positional window of validity,
 *     requiring splitting into multiple validity intervals
 *
 * If we collapse a branch as soon as pos >= first_kmer_end_position of root then:
 * - n in G(pos) since last_kmer_start_position(n) < last_kmer_end_position(n) < first_kmer_start_position(n) <= pos
 * - break n length-wise
 *   - need to prove that the partially constructed nodes adjacent to n will reference the correct split
 *    - yes: neighbours of split before n' are fully defined 
 *    - yes: neighbours of n' are fully defined
 *    - yes: neighbours of post n' are partially defined iff they only connect to post n' -> connect to post n' post-split
 *  
 * Leaf collapse:
 *    * - - - - -   <- leaf     
 *                \    
 * * - * - * - * - * - ?
 *                 ^
 *                root
 * 
 * Conditions for backward leaf collapse match that of path collapse when the
 * root node is considered to be the node adjacent to the leaf
 * 
 *     - - - - - *      <- collapse this leave into a main path
 *   /
 * * - * - * - * - * - ?
 * 
 * Similarly, forward leaf collapse requires fully defined path nodes up to last_kmer_end_position(leaf)
 * 
 * For simplicity, each node processed for all collapse types once:
 *                             |----maxSupportWidth+maxPathLength----|
 *    |----maxcollapseLength---|                                     |----maxcollapseLength---|
 *    ^                                                              ^                        ^
 *    |                                      ^^^^                 process                     |
 * unchanged                        node being processed          offset                  input position 
 *  offset
 * 
 * Note: the output graph is not minimal and may contain adjacent nodes that could be merged
 * 
 * @author Daniel Cameron
 *
 */
public class PathCollapseIterator extends CollapseIterator implements DeBruijnGraph<KmerPathSubnode> {
	private final boolean bubblesAndLeavesOnly;
	private final double minimumPathNodeEntropy;
	public PathCollapseIterator(
			Iterator<KmerPathNode> it,
			int k,
			int maxPathCollapseLength,
			int maxBasesMismatch,
			boolean bubblesAndLeavesOnly,
			double minimumPathNodeEntropy) {
		super(it, k, maxPathCollapseLength, maxBasesMismatch, 0, 0);
		this.bubblesAndLeavesOnly = bubblesAndLeavesOnly;
		this.minimumPathNodeEntropy = minimumPathNodeEntropy;
	}
	@Override
	protected boolean reprocessMergedNodes() {
		// no need to reprocess since we exhaustively check paths
		return false;
	}
	protected boolean collapse(KmerPathNode node, int maxCollapseLength) {
		KmerPathSubnode root = new KmerPathSubnode(node);
		List<KmerPathSubnode> nextNodes = root.next();
		for (int i = 0; i < nextNodes.size(); i++) {
			for (int j = i + 1; j < nextNodes.size(); j++) {
				if (collapseSimilarPath(node, nextNodes.get(i), nextNodes.get(j), true, true, true, maxCollapseLength)) return true;
			}
		}
		List<KmerPathSubnode> prevNodes = root.prev();
		for (int i = 0; i < prevNodes.size(); i++) {
			for (int j = i + 1; j < prevNodes.size(); j++) {
				if (collapseSimilarPath(node, prevNodes.get(i), prevNodes.get(j), true, false, false, maxCollapseLength)) return true;
			}
		}
		return false;
	}
	private boolean collapseSimilarPath(KmerPathNode root, KmerPathSubnode startNodeA, KmerPathSubnode startNodeB, boolean findLeaf, boolean findCommonChild, boolean traverseForward, int maxCollapseLength) {
		if (!hasSufficientEntropy(startNodeA.node())) return false;
		if (!hasSufficientEntropy(startNodeB.node())) return false;
		KmerPathNodePath pathA = new KmerPathNodePath(startNodeA, traverseForward, maxCollapseLength);
		KmerPathNodePath pathB = new KmerPathNodePath(startNodeB, traverseForward, maxCollapseLength);
		if (pathA.pathLength() <= maxCollapseLength && pathB.pathLength() <= maxCollapseLength) {
			Set<KmerPathNode> set = Collections.newSetFromMap(new IdentityHashMap<KmerPathNode, Boolean>());
			set.add(root);
			set.add(pathA.headPath());
			set.add(pathB.headPath());
			if (set.size() == 3) {
				return collapseSimilarPath(set, pathA, pathB, pathBasesDifferent(pathA, pathB, traverseForward), findLeaf, findCommonChild, traverseForward);
			}
		}
		return false;
	}
	/**
	 * Recusive node traversal of path trees looking for similar paths
	 * 
	 * Need to simultaneous traverse across both trees comparing all possible path combinations until a match is found. 
	 * 
	 * Invariant: pathA and pathB are unchanged when returning false
	 * 
	 * Although the underlying graph is a DAG, KmerPathSubnode traversal is not. The same
	 * underlying KmerPathNode can be present multiple times in either path. Example scenario:
	 * [1,10] AAAA -> [2,11] AAAA -> [3,12] AAAG = [1,10] AAAAG 
	 * [1,10] AAAA -> [2,11] AAAG -> [3,12] AAGG = [1,10] AAAGG
	 * only 1 kmer difference so should collapse but when merging, we are unable to fragment
	 * [1,11] AAAA in such a way that each KmerPathSubnode contains a single KmerPathNode
	 * per KmerPathSubnode, we have to fragment the KmerPathNode across all boundaries.
	 * For now, we handle this by not collapsing if we encounter a KmerPathNode repeat. 
	 * 
	 */
	private boolean collapseSimilarPath(
			final Set<KmerPathNode> onPath,
			final KmerPathNodePath pathA,
			final KmerPathNodePath pathB,
			final int basesDifferent,
			final boolean findLeaf,
			final boolean findCommonChild,
			final boolean traverseForward) {
		nodesTraversed++;
		if (!couldMatch(pathA, pathB, basesDifferent, traverseForward)) return false;
		if (findLeaf) {
			if (tryLeafCollapse(pathA, pathB, traverseForward)) return true;
			if (tryLeafCollapse(pathB, pathA, traverseForward)) return true; 
		}
		int pathAlength = pathA.pathLength();
		int pathBlength = pathB.pathLength();
		if (pathA.pathLength() <= pathB.pathLength()) {
			while (pathA.dfsNextChild()) {
				KmerPathNode added = pathA.headPath();
				int additionalBasesDifferent = headNodeBasesDifferent(pathA, pathB, traverseForward);
				if (!onPath.contains(added) && hasSufficientEntropy(added)) {
					onPath.add(added);
					pathB.dfsResetChildTraversal();
					if (collapseSimilarPath(onPath, pathA, pathB, basesDifferent + additionalBasesDifferent, findLeaf, findCommonChild, traverseForward)) return true;
					assert(added == pathA.headPath());
					onPath.remove(added);
				} else if (added == pathB.headPath() && findCommonChild && couldMatch(pathA, pathB, basesDifferent + additionalBasesDifferent, traverseForward)) {
					if (tryPathCollapse(pathA, pathB, traverseForward)) return true;
				}
				pathA.pop();
				assert(pathAlength == pathA.pathLength());
				assert(pathBlength == pathB.pathLength());
			}
			//pathA.dfsPop(); // done with this node
		} else {
			while (pathB.dfsNextChild()) {
				KmerPathNode added = pathB.headPath();
				int additionalBasesDifferent = headNodeBasesDifferent(pathB, pathA, traverseForward);
				if (!onPath.contains(added) && hasSufficientEntropy(added)) {
					onPath.add(added);
					if (collapseSimilarPath(onPath, pathA, pathB, basesDifferent + additionalBasesDifferent, findLeaf, findCommonChild, traverseForward)) return true;
					assert(added == pathB.headPath());
					onPath.remove(added);
				} else if (added == pathA.headPath() && findCommonChild && couldMatch(pathA, pathB, basesDifferent + additionalBasesDifferent, traverseForward)) {
					if (tryPathCollapse(pathA, pathB, traverseForward)) return true;
				}
				pathB.pop();
				assert(pathAlength == pathA.pathLength());
				assert(pathBlength == pathB.pathLength());
			}
		}
		assert(pathA.pathLength() == pathAlength);
		assert(pathB.pathLength() == pathBlength);
		return false;
	}
	private boolean couldMatch(KmerPathNodePath pathA, KmerPathNodePath pathB, int basesDifferent, boolean traverseForward) {
		if (Defaults.SANITY_CHECK_ASSEMBLY_GRAPH) {
			// can't actually do this assertion as the path difference may actually be violated when a non-bubble path is merged
			//      F - -
			//     /     \
			//    D - E   \
			//   /     \   \
			//  A - B - C - G
			//
			// Step 1) DE merged into B
			//      F - -
			//     /     \
			//   B1'-B2   \
			//   /     \   \
			//  A       C - G
			//
			// 
			// Transition from B1' to F is now an invalid kmer path
			// number of bases different now depends on which base call is 
			assert(!bubblesAndLeavesOnly || basesDifferent == pathBasesDifferent(pathA, pathB, traverseForward));
		}
		if (basesDifferent > maxBasesMismatch) return false;
		// paths don't share any common interval - no way for them to be the same length and still overlap
		if (!pathsOverlap(pathA, pathB)) return false;
		return true;
	}
	private int pathBasesDifferent(KmerPathNodePath pathA, KmerPathNodePath pathB, boolean traverseForward) {
		int basesDifference;
		if (traverseForward) {
			basesDifference = DeBruijnSequenceGraphNodeUtil.basesDifferent(k, pathA.currentPath(), pathB.currentPath()); 
		} else {
			basesDifference = DeBruijnSequenceGraphNodeUtil.reverseBasesDifferent(k, pathA.currentPath(), pathB.currentPath());
		}
		return basesDifference;
	}
	/**
	 * Number of bases that the headPath head node differs from the refPath
	 */
	private int headNodeBasesDifferent(final KmerPathNodePath headPath, final KmerPathNodePath refPath, final boolean traverseForward) {
		int diff = 0;
		final TraversalNode headNode = headPath.headNode();
		final KmerPathNode headPathNode = headPath.headPath();
		assert(headNode.parent() != null);
		final int headPathLengthExcludingHead = headNode.parent().pathLength();
		final int headPathStartOffset = headNode.pathLength() - headPathNode.length();
		final int headPathEndOffset = headNode.pathLength() - 1;
		TraversalNode refNode = refPath.headNode();
		while (refNode != null && refPath.pathLength() > headPathLengthExcludingHead) {
			final KmerPathNode refPathNode = refNode.node().node();
			final int nodePathStartOffset = refNode.pathLength() - refPathNode.length();
			final int nodePathEndOffset = refNode.pathLength() - 1;
			final int commonPathStartOffset = Math.max(headPathStartOffset, nodePathStartOffset);
			final int commonPathEndOffset  = Math.min(headPathEndOffset, nodePathEndOffset);
			if (commonPathEndOffset < commonPathStartOffset) break;
			for (int offset = commonPathStartOffset; offset <= commonPathEndOffset; offset++) {
				if ((traverseForward && !KmerEncodingHelper.lastBaseMatches(k,
						headPathNode.kmer(offset - headPathStartOffset),
						refPathNode.kmer(offset - nodePathStartOffset)))
					|| (!traverseForward && !KmerEncodingHelper.firstBaseMatches(k,
						headPathNode.kmer(headPathNode.length() - 1 - (offset - headPathStartOffset)),
						refPathNode.kmer(refPathNode.length() - 1 - (offset - nodePathStartOffset))))) {
					diff++;
				}
			}
			refNode = refNode.parent();
		}
		return diff;
	}
	private boolean pathsOverlap(KmerPathNodePath pathA, KmerPathNodePath pathB) {
		boolean overlaps = IntervalUtil.overlapsClosed(
				pathA.headNode().startPositionOfAnchorKmer(), pathA.headNode().endPositionOfAnchorKmer(),
				pathB.headNode().startPositionOfAnchorKmer(), pathB.headNode().endPositionOfAnchorKmer());
		return overlaps;
	}
	private boolean tryPathCollapse(KmerPathNodePath pathA, KmerPathNodePath pathB, boolean traverseForward) {
		if (pathA.headPath() == pathB.headPath()) {
			if (pathA.pathLength() == pathB.pathLength()) {
				// remove common trailing node
				pathA.pop();
				pathB.pop();
				assert(repeatedKmerPathNodeCount(pathA, pathB) == 0);
				assert(pathsOverlap(pathA, pathB));
				List<KmerPathSubnode> lA = new ArrayList<KmerPathSubnode>(pathA.headNode().overlapping(pathB.headNode()).asSubnodes());
				List<KmerPathSubnode> lB = new ArrayList<KmerPathSubnode>(pathB.headNode().overlapping(pathA.headNode()).asSubnodes());
				if (pathA.pathWeight() < pathB.pathWeight()) {
					if (!bubblesAndLeavesOnly || isBubblePath(lA)) {
						merge(lA, lB, 0, 0);
						branchesCollapsed++;
						return true;
					}
				} else {
					if (!bubblesAndLeavesOnly || isBubblePath(lB)) {
						merge(lB, lA, 0, 0);
						branchesCollapsed++;
						return true;
					}
				}
			}
		}
		return false;
	}
	/**
	 * A path is considered a bubble if it each node only has a single source and successor 
	 * @param path path excluding starting node of bubble, but including ending node
	 * @return true if the path is a bubble, false otherwise
	 */
	private boolean isBubblePath(List<KmerPathSubnode> path) {
		for (int i = 0; i < path.size() - 1; i++) {
			KmerPathSubnode sn = path.get(i);
			if (sn.next().size() != 1) return false;
			if (sn.prev().size() != 1) return false;
		}
		return true;
	}
	private int repeatedKmerPathNodeCount(KmerPathNodePath... paths) {
		Set<KmerPathNode> set = Collections.newSetFromMap(new IdentityHashMap<KmerPathNode, Boolean>());
		int nodeCount = 0;
		for (KmerPathNodePath path : paths) {
			nodeCount += path.currentPath().size();
			set.addAll(path.currentPath());
		}
		// we if shrunk in size then we have a repeat
		return nodeCount - set.size();
	}
	private boolean tryLeafCollapse(KmerPathNodePath leaf, KmerPathNodePath path, boolean traverseForward) {
		// leaf can't be longer than the path
		if (leaf.pathLength() > path.pathLength()) return false;
		if (leaf.pathWeight() > path.pathWeight()) return false;
		TraversalNode firstLeaf = leaf.headNode().overlapping(path.headNode()).firstTerminalLeaf();
		if (firstLeaf == null) return false;
		assert(repeatedKmerPathNodeCount(leaf, path) == 0);
		assert(pathsOverlap(leaf, path));
		int leafSkip = 0;
		int pathSkip = traverseForward ? 0 : path.pathLength() - leaf.pathLength();
		merge(new ArrayList<KmerPathSubnode>(firstLeaf.asSubnodes()),
				new ArrayList<KmerPathSubnode>(path.headNode().overlapping(firstLeaf).asSubnodes()),
				leafSkip, pathSkip);
		leavesCollapsed++;
		return true;
	}
	private boolean hasSufficientEntropy(KmerPathNode node) {
		if (minimumPathNodeEntropy <= 0) return true;
		double entropy = SequenceUtil.shannonEntropy(KmerEncodingHelper.baseCounts(k, node.pathKmers()));
		return entropy > minimumPathNodeEntropy;
	}
	@Override
	public int getWeight(KmerPathSubnode node) {
		return node.weight();
	}
	@Override
	public List<KmerPathSubnode> next(KmerPathSubnode node) {
		return node.next();
	}
	@Override
	public List<KmerPathSubnode> prev(KmerPathSubnode node) {
		return node.prev();
	}
	@Override
	public void removeNode(KmerPathSubnode node) {
		throw new UnsupportedOperationException();
	}
	@Override
	public void removeEdge(KmerPathSubnode source, KmerPathSubnode sink) {
		throw new UnsupportedOperationException();
	}
	@Override
	public void addNode(KmerPathSubnode node) {
		throw new UnsupportedOperationException();
	}
	@Override
	public void addEdge(KmerPathSubnode source, KmerPathSubnode sink) {
		throw new UnsupportedOperationException();
	}
	@Override
	public Collection<KmerPathSubnode> allNodes() {
		throw new UnsupportedOperationException();
	}
	@Override
	public String toString(Iterable<? extends KmerPathSubnode> path) {
		return String.format("[%d-%d] %s",
			path.iterator().next().firstStart(),
			path.iterator().next().firstEnd(),
			new String(KmerEncodingHelper.baseCalls(DeBruijnSequenceGraphNodeUtil.asKmers(path), k)));
	}
	@Override
	public int getK() {
		return k;
	}
	@Override
	public long getKmer(KmerPathSubnode node) {
		return node.node().lastKmer();
	}
	@Override
	public boolean isReference(KmerPathSubnode node) {
		return node.node().isReference();
	}
}