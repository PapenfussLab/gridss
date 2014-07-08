package au.edu.wehi.idsv.debruijn;

import htsjdk.samtools.util.Log;

import java.util.ArrayDeque;
import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Set;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Ordering;
import com.google.common.collect.Sets;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;

/**
 * Compressed De Bruijn graph in which each node is a kmer path with no branches 
 * @author Daniel Cameron
 *
 */
public abstract class DeBruijnPathGraph<T extends DeBruijnNodeBase, PN extends PathNode<T>> {
	protected abstract PN createPathNode(LinkedList<Long> path, DeBruijnGraphBase<T> graph);
	protected DeBruijnGraphBase<T> graph;
	protected Set<PN> paths;
	protected Map<PN, List<PN>> pathNext;
	protected Map<PN, List<PN>> pathPrev;
	private static Log log = Log.getInstance(DeBruijnPathGraph.class);
	public DeBruijnPathGraph(DeBruijnGraphBase<T> graph, long seed) {
		this.graph = graph;
		generatePathGraph(seed);
	}
	public Set<PN> getPaths() { return paths; }
	/**
	 * generates a path graph from the given seed kmer
	 * @param seed
	 */
	private void generatePathGraph(long seed) {
		Queue<Long> frontier = new ArrayDeque<Long>();
		frontier.add(seed);
		Map<Long, PN> pathStart = Maps.newHashMap();
		Map<Long, PN> pathEnd = Maps.newHashMap();
		paths = Sets.newHashSet();
		while (!frontier.isEmpty()) {
			long kmer = frontier.poll();
			if (pathStart.containsKey(kmer)) continue;
			if (pathEnd.containsKey(kmer)) continue;
			PN path = traverseBranchless(kmer);
			for (long adj : graph.prevStates(path.getFirst(), null, null)) {
				frontier.add(adj);
			}
			for (long adj : graph.nextStates(path.getLast(), null, null)) {
				frontier.add(adj);
			}
			// Add path to graph
			pathStart.put(path.getFirst(), path);
			pathEnd.put(path.getLast(), path);
			paths.add(path);
		}
		pathNext = Maps.newHashMap();
		pathPrev = Maps.newHashMap();
		for (PN path : paths) {
			// construct edges
			List<Long> prevKmers = graph.prevStates(path.getFirst(), null, null);
			List<PN> prev = Lists.newArrayListWithExpectedSize(prevKmers.size());
			for (long adj : prevKmers) {
				prev.add(pathEnd.get(adj));
			}
			List<Long> nextKmers = graph.nextStates(path.getLast(), null, null);
			List<PN> next = Lists.newArrayListWithExpectedSize(nextKmers.size());
			for (long adj : nextKmers) {
				next.add(pathStart.get(adj));
			}
		}
	}
	/**
	 * Returns all paths that follow the given path
	 * @param path path
	 * @return successor paths
	 */
	public List<PN> nextPath(PathNode path) {
		return pathNext.get(path);
	}
	/**
	 * Returns all paths that follow the given path
	 * @param path path
	 * @return preceeding paths
	 */
	public List<PN> prevPath(PathNode path) {
		return pathPrev.get(path);
	}
	/**
	 * Self-intersecting paths will never be traversed
	 * by a contig generated from simple kmer path
	 * so can be removed from the graph
	 */
	public void removeSelfIntersectingPaths() {
		for (PN path : paths) {
			if (nextPath(path).contains(path)) {
				nextPath(path).remove(path);
			}
			if (prevPath(path).contains(path)) {
				prevPath(path).remove(path);
			}
		}
	}
	/**
	 * Collapses all paths that differ by at most maxDifference bases
	 * into the path with the highest weight.
	 * @param maxDifference maximum differences along the path. Note that by collapsing
	 * a graph with n differences, alternate paths that share a PathNode  may have up to n differences
	 * from their original due to a differences merging the paths
	 *       E      F                    E             F
	 *        \    /                      \           / 
	 *    G - H - I - J    ->              \         / 
	 *   /              \                   \       /        
	 *  A  -  B  -  C -  D          A - B1 - B2 - C1 - C2 - D
	 *  
	 *  E-H-I-J is turned in E-B1-B2-C1-C2-F resulting in a different path
	 *  @param bubblesOnly only collapse bubbles. No other paths will be affected
	 */
	public void collapseSimilarPaths(int maxDifference, boolean bubblesOnly) {
		// Keep collapsing until we can't anymore
		boolean collapsed = true;
		while (collapsed) {
			collapsed = false;
			for (PathNode start : paths) {
				// TODO: collapse the path with the weakest support first
				while (collapsePaths(maxDifference, bubblesOnly, start)) {
					collapsed = true;
				}
			}
		}
	}
	private boolean collapsePaths(int maxDifference, boolean bubblesOnly, PathNode start) {
		LinkedList<PN> listA = Lists.newLinkedList();
		LinkedList<PN> listB = Lists.newLinkedList();
		List<PN> next = nextPath(start);
		for (int i = 0; i < next.size() - 1; i++) {
			for (int j = i + 1; j < next.size(); j++) {
				listA.add(next.get(i));
				listB.add(next.get(j));
				if (collapsePaths(maxDifference, bubblesOnly, listA, listB, next.get(i).length(), next.get(j).length())) {
					return true;
				}
				assert(listA.size() == 1);
				assert(listB.size() == 1);
				listA.removeFirst();
				listB.removeFirst();
			}
		}
		return false;
	}
	private boolean collapsePaths(
			int differencesAllowed,
			boolean bubblesOnly,
			LinkedList<PN> pathA,
			LinkedList<PN> pathB,
			int pathALength,
			int pathBLength) {
		// paths have diverged too far
		if (basesDifferent(pathA, pathB) > differencesAllowed) return false;
		
		if (shareNextPaths(pathA.getLast(), pathB.getLast())) {
			// only accept base difference, not indels
			if (pathALength != pathBLength) return false;
			
			// we have a common path!
			PN common = pathA.getLast();
			pathA.removeLast();
			pathB.removeLast();
			if (bubblesOnly && (!isBubble(pathA) && !isBubble(pathB))) return false;
			boolean collapsed = mergePath(pathA, pathB);
			pathA.addLast(common);
			pathB.addLast(common);
			return collapsed;
		}
		// short circuit if only processing bubbles
		if (bubblesOnly && (!isBubble(pathA) && !isBubble(pathB))) return false;
		
		if (pathALength <= pathBLength) {
			for (PN next : nextPath(pathA.getLast())) {
				pathA.addLast(next);
				if (collapsePaths(differencesAllowed, bubblesOnly, pathA, pathB, pathALength + next.length(), pathBLength)) {
					return true;
				}
				pathA.removeLast();
			}
		} else {
			for (PN next : nextPath(pathB.getLast())) {
				pathB.addLast(next);
				if (collapsePaths(differencesAllowed, bubblesOnly, pathA, pathB, pathALength, pathBLength + next.length())) {
					return true;
				}
				pathB.removeLast();
			}
		}
		return false; 
	}
	private boolean shareNextPaths(PN pathA, PN pathB) {
		// Can't just check if the end up at the same kmer as one of the
		// paths may have been collasped into another kmer path 
		// return KmerEncodingHelper.nextState(graph.getK(), pathA.getLast(), (byte)'A') == KmerEncodingHelper.nextState(graph.getK(), pathA.getLast(), (byte)'A');
		// Check the next states are the same
		// Since next states are ordered, we can just compare that the two lists are the same 
		return Iterables.elementsEqual(nextPath(pathA), nextPath(pathB));
	}
	private int basesDifferent(LinkedList<PN> pathA, LinkedList<PN> pathB) {
		Iterator<Long> itA = PathNode.kmerIterator(pathA);
		Iterator<Long> itB = PathNode.kmerIterator(pathB);
		int differences = KmerEncodingHelper.basesDifference(graph.getK(), itA.next(), itB.next());
		while (itA.hasNext() && itB.hasNext()) {
			if (KmerEncodingHelper.lastBaseMatches(graph.getK(), itA.next(), itB.next())) {
				differences++;
			}
		}
		return differences;
	}
	/**
	 * Merges the given paths together
	 * @param pathA first path to merge 
	 * @param pathB second path to merge
	 * @return true if a merge could be performed, false otherwise
	 */
	public boolean mergePath(LinkedList<PN> pathA, LinkedList<PN> pathB) {
		int weightA = getWeight(pathA);
		int weightB = getWeight(pathB);
		boolean isBubbleA = isBubble(pathA);
		boolean isBubbleB = isBubble(pathB);
		if (isBubbleA && isBubbleB) {
			if (weightA < weightB) {
				return mergeBubble(pathA.getFirst(), pathB);
			} else {
				return mergeBubble(pathB.getFirst(), pathA);
			}
		} else if (isBubbleA) {
			return mergeBubble(pathA.getFirst(), pathB);
		} else if (isBubbleB) {
			return mergeBubble(pathA.getFirst(), pathB);
		}
		if (bigBubbleWarned < 3) {
			log.warn(String.format("Encountered path requiring merge but could not be merged. (This message will not be repeated)"));
			pathMergingWarned++;
		}
		return false;
	}
	private static int pathMergingWarned = 0;
	/**
	 * Merges the given paths together
	 * @param pathA first path to merge 
	 * @param pathB second path to merge
	 */
	public boolean mergeBubble(PN bubble, LinkedList<PN> path) {
		int bubbleWeight = bubble.getWeight();
		int pathWeight = getWeight(path);
		if (bubbleWeight > pathWeight) {
			if (bigBubbleWarned < 3) {
				log.warn(String.format("Failed to collapse bubble: weight of kmer bubble is greater than the path weight. (This message will not be repeated)"));
				bigBubbleWarned++;
			}
			return false;
		}
		// replace edges that point to the bubble to point to the path
		for (PN prev : prevPath(bubble)) {
			listReplace(nextPath(prev), bubble, path.getFirst());
		}
		for (PN next : nextPath(bubble)) {
			listReplace(prevPath(next), bubble, path.getLast());
		}
		// Remove the bubble path
		this.pathNext.remove(bubble);
		this.pathPrev.remove(bubble);
		this.paths.remove(bubble);
		
		// Update kmers in the underlying graph
		Iterator<Long> itBubble = bubble.getPath().iterator();
		for (PathNode p : path) {
			LinkedList<Long> altPath = Lists.newLinkedList();
			for (int i = 0; i < p.length(); i++) {
				if (!itBubble.hasNext()) throw new RuntimeException("Sanity check failure: bubble is shorter than path being merged to.");
				altPath.add(itBubble.next());
			}
			p.merge(altPath, graph);
		}
		return true;
	}
	private void listReplace(List<PN> list, PN toReplace, PN replaceWith) {
		if (list.contains(toReplace)) {
			if (list.contains(replaceWith)) {
				// contains both: just remove
				list.remove(toReplace);
			} else {
				// only contains toReplace, just overwrite with the new value
				list.set(list.indexOf(toReplace), replaceWith);
			}
		}
	}
	private static int bigBubbleWarned = 0;
	/**
	 * Gets the total weight of all kmers on the given path
	 * @param path path
	 * @return total weight
	 */
	public int getWeight(LinkedList<PN> path) {
		int weight = 0;
		for (PathNode n : path) {
			weight += n.getWeight();
		}
		return weight;
	}
	/**
	 * A bubble is a de bruijn graph path diverges from a reference path
	 * (usually due to a sequencing error), then converges back to the reference
	 * 
	 * 
	 * @param path path to test
	 * @return true if the path could be a bubble, false otherwise
	 */
	public boolean isBubble(LinkedList<PN> path) {
		return path.size() == 1
				&& prevPath(path.getFirst()).size() == 1 
				&& nextPath(path.getLast()).size() == 1;
	}
	/**
	 * Traverses kmers until a branch is found
	 * @param seed starting kmer
	 * @return unique path
	 */
	protected PN traverseBranchless(long seed) {
		LinkedList<Long> path = new LinkedList<Long>();
		Set<Long> visited = Sets.newHashSet();
		path.add(seed);
		for(List<Long> adj = graph.nextStates(path.getLast(), null, null); adj.size() == 1 && graph.prevStates(adj.get(0), null, null).size() <= 1; adj = graph.nextStates(path.getLast(), null, null)) {
			if (visited.contains(adj)) {
				// circular contig
			}
			path.addLast(adj.get(0));
			visited.add(adj.get(0));
		}
		for(List<Long> adj = graph.prevStates(path.getFirst(), null, null); adj.size() == 1 && graph.nextStates(adj.get(0), null, null).size() <= 1; adj = graph.prevStates(path.getFirst(), null, null)) {
			if (visited.contains(adj)) {
				// circular contig
				break;
			}
			path.addFirst(adj.get(0));
			visited.add(adj.get(0));
		}
		return createPathNode(path, graph);
	}
	/**
	 *  Greedily traverse a path based on the given choice ordering
	 * @param startNode starting position
	 * @param forwardChoice ordering when traversing forward
	 * @param backwardChoice ordering when traversing backward 
	 * @return traversal
	 */
	public LinkedList<PN> greedyTraverse(PN startNode, Comparator<PN> forwardChoice, Comparator<PN> backwardChoice) {
		LinkedList<PN> path = new LinkedList<PN>();
		Set<PN> visited = Sets.newHashSet();
		path.add(startNode);
		visited.add(startNode);
		// assemble back
		PriorityQueue<PN> prevCandidates = new PriorityQueue<PN>(backwardChoice);
		for (List<PN> nodeList = prevPath(path.getFirst()); ; nodeList = prevPath(path.getFirst())) {
			prevCandidates.clear();
			for (PN node : nodeList) {
				if (visited.contains(node)) continue; // no loops
				prevCandidates.add(node);
			}
			// we're done
			if (prevCandidates.isEmpty()) break;
			path.addFirst(prevCandidates.poll());
		}
		// assemble forward
		PriorityQueue<PN> nextCandidates = new PriorityQueue<PN>(forwardChoice);
		for (List<PN> nodeList = nextPath(path.getLast()); ; nodeList = nextPath(path.getLast())) {
			nextCandidates.clear();
			for (PN node : nodeList) {
				if (visited.contains(node)) continue; // no loops
				nextCandidates.add(node);
			}
			// we're done
			if (nextCandidates.isEmpty()) break;
			path.addLast(nextCandidates.poll());
		}
		return path;
	}
	public PN getNodeContaining(long kmer) {
		for (PN node : paths) {
			if (node.contains(kmer)) return node;
		}
		return null;
	}
	/**
	 * Ordering of kmers by kmer weight. 
	 */
	// Needs to be non-static to access graph.getK()
	public Ordering<PN> ByAverageKmerWeightDesc = new Ordering<PN>() {
		public int compare(PN o1, PN o2) {
			return Doubles.compare(o1.getWeight() / (o1.length() + graph.getK() - 1), o2.getWeight() / (o2.length() + graph.getK() - 1));
		}
	}.reverse();
	/**
	 * Ordering of kmers by kmer weight. 
	 */
	// Needs to be non-static as we can't have a generic type static property
	public Ordering<PN> ByPathTotalWeightDesc = new Ordering<PN>() {
		public int compare(PN o1, PN o2) {
			return Ints.compare(o1.getWeight(), o2.getWeight());
		}
	}.reverse();
	/**
	 * Ordering of kmers by kmer weight. 
	 */
	public Ordering<PN> ByMaxKmerWeightDesc = new Ordering<PN>() {
		public int compare(PN o1, PN o2) {
			return Ints.compare(o1.getMaxKmerWeight(), o2.getMaxKmerWeight());
		}
	}.reverse();
}
