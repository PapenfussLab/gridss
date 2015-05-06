package au.edu.wehi.idsv.graph;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import au.edu.wehi.idsv.util.FirstOverflowList;

import com.google.common.base.Function;
import com.google.common.collect.Iterators;

public class PathNode<T> {
	/**
	 * Sequential node identifier within this path graph.
	 * Used improve lookup performance
	 * @return node identifier
	 */
	public int getNodeId() {
		return nodeId;
	}
	/**
	 * Performance optimisation: offset of this PathNode within the containing PathGraph lookups
	 */
	void setNodeId(int nodeId) {
		this.nodeId = nodeId;
	}
	private int nodeId = -1;
	private int weight = 0;
	private ArrayList<T> path;
	private ArrayList<LinkedList<T>> equivalents = null;
	private FirstOverflowList<T> fol; // TODO: refactor and inherit FirstOverflowList
	/**
	 * Creates a new path graph node from the given node sequence
	 * @param path node sequence
	 * @param graph parent graph
	 */
	public PathNode(Collection<T> nodes, WeightedDirectedGraph<T> graph) {
		assert(nodes != null);
		assert(graph != null);
		assert(nodes.size() > 0);
		this.path = new ArrayList<T>(nodes.size());
		this.path.addAll(nodes);
		for (T n : path) {
			this.weight += graph.getWeight(n);
		}
		this.fol = new FirstOverflowList<T>(path, equivalents);
	}
	/**
	 * Merges the given path into the current path node
	 * @param path node sequence
	 * @param graph parent graph
	 * @return 
	 */
	public static <T> PathNode<T> concat(Iterable<? extends PathNode<T>> path, WeightedDirectedGraph<T> graph) {
		int length = nodeLength(path);
		return new PathNode<T>(path, 0, length, graph);
	}
	/**
	 * Concatenates the given path nodes together
	 * @param nodes path
	 * @param startOffset number of nodes to skip
	 * @param length length of path to take
	 * @param graph graph
	 */
	public PathNode(Iterable<? extends PathNode<T>> nodes, int startOffset, int length, WeightedDirectedGraph<T> graph) {
		path = new ArrayList<T>(length);
		for (PathNode<T> pn : nodes) {
			if (pn.equivalents != null) {
				equivalents = new ArrayList<LinkedList<T>>(length);
				fol = new FirstOverflowList<T>(path, equivalents);
				break;
			}
		}
		int offset = 0;
		for (PathNode<T> pn : nodes) {
			for (int j = 0; j < pn.length() && offset - startOffset < length; j++, offset++) {
				if (offset < startOffset) {
					// skip
				} else {
					T n = pn.path.get(j);
					path.add(n);
					if (equivalents != null) {
						if (pn.equivalents != null && pn.equivalents.get(j) != null) {
							equivalents.add(new LinkedList<T>(pn.equivalents.get(j)));
						} else {
							equivalents.add(null);
						}
					}
				}
			}
		}
		fol = new FirstOverflowList<T>(path, equivalents);
		recalculateWeight(graph);
	}
	private void recalculateWeight(WeightedDirectedGraph<T> graph) {
		weight = 0;
		for (T n : path) {
			weight += graph.getWeight(n);
		}
		if (equivalents != null) {
			for (LinkedList<T> eq : equivalents) {
				if (eq != null) {
					for (T n : eq) {
						weight += graph.getWeight(n);
					}
				}
			}
		}
	}
	public List<T> getPath() { return path; }
	/**
	 * All nodes on this path including nodes that have been merged from similar paths
	 * @return iterator over all nodes
	 */
	public Iterable<T> allNodes() {
		return fol.asFlattenedIterable();
	}
	public List<List<T>> getPathAllNodes() {
		return fol;
	}
	/**
	 * Gets node paths that have been merged into this path
	 * @return
	 */
	public T first() { return path.get(0); }
	public T last() { return path.get(length() - 1); }
	public int length() { return path.size(); }
	public int nodes() {
		int nodeCount = path.size();
		if (equivalents != null) {
			for (LinkedList<T> l : equivalents) {
				if (l != null) {
					nodeCount += l.size();
				}
			}
		}
		return nodeCount;
	}
	public int weight() {
		return weight;
	}
	private int weight(int offset, WeightedDirectedGraph<T> graph) {
		int w = graph.getWeight(path.get(offset));
		if (equivalents != null) {
			LinkedList<T> c = equivalents.get(offset);
			if (c != null) {
				for (T n : c) {
					w += graph.getWeight(n);
				}
			}
		}
		return w;
	}
	public int maxNodeWeight(WeightedDirectedGraph<T> graph) {
		int maxNodeWeight = 0;
		for (Collection<T> c : getPathAllNodes()) {
			int nodeWeight = 0;
			for (T n : c) {
				nodeWeight += graph.getWeight(n);
			}
			maxNodeWeight = Math.max(maxNodeWeight, nodeWeight); 
		}
		return maxNodeWeight;
	}
	/**
	 * Merges the nodes on this alternate path into this path
	 * @param alternatePath alternate path to merge
	 */
	public void merge(Iterable<? extends PathNode<T>> alternatePath, WeightedDirectedGraph<T> graph) {
		assert(nodeLength(alternatePath) == length());
		merge(alternatePath, 0, this.length(), 0, graph);
	}
	/**
	 * Merges the nodes on this alternate path into this path
	 * @param alternatePath alternate path to merge
	 * @param alternatePathKmerOffset number of initial nodes on alternate path to skip before merging
	 * @param alternatePathKmers number of nodes on alternate path to merge
	 * @param nodeOffset position to start merging nodes
	 * @param length alternatePath alternate path to merge
	 */
	public void merge(
			Iterable<? extends PathNode<T>> alternatePath,
			int alternatePathKmerOffset,
			int alternatePathKmers,
			int nodeOffset,
			WeightedDirectedGraph<T> graph) {
		assert(nodeLength(alternatePath) >= alternatePathKmers + alternatePathKmerOffset);
		if (alternatePathKmers > length() - nodeOffset) throw new IllegalArgumentException("Alternate path is too long.");
		for (PathNode<T> pn : alternatePath) {
			for (int j = 0; j < pn.length() && alternatePathKmers > 0; j++) {
				if (alternatePathKmerOffset > 0) {
					alternatePathKmerOffset--;
				} else {
					merge(nodeOffset++, pn, j, graph);
					alternatePathKmers--;
				}
			}
		}
	}
	/**
	 * Merges the nodes on this alternate path into this path
	 * @param alternatePath alternate path to merge
	 * @param offset number of initial nodes that the alternate path does not contain.
	 * The first alternate path node will be merged into the node at offset position.
	 * @param length 
	 
	public void merge(Iterable<? extends PathNode<T>> alternatePath, int offset, WeightedDirectedGraph<T> graph) {
		if (nodeLength(alternatePath) + offset > length()) throw new IllegalArgumentException("Alternate path is too long.");
		int i = 0;
		for (PathNode<T> node : alternatePath) {
			for (int j = 0; j < node.length(); j++) {
				merge(offset + i++, node, j, graph);
			}
		}
	}*/
	protected void merge(int offset, PathNode<T> pn, int pnOffset, WeightedDirectedGraph<T> graph) {
		assert(offset >= 0 && offset < length());
		assert(pnOffset >= 0 && pnOffset < pn.length());
		if (equivalents == null) {
			equivalents = new ArrayList<LinkedList<T>>(length());
			fol = new FirstOverflowList<T>(path, equivalents);
			for (int i = 0; i < length(); i++) {
				equivalents.add(null);
			}
		}
		LinkedList<T> eq = equivalents.get(offset);
		if (eq == null) {
			eq = new LinkedList<T>();
			equivalents.set(offset, eq);
		}
		eq.add(pn.path.get(pnOffset));
		weight += graph.getWeight(pn.path.get(pnOffset));
		if (pn.equivalents != null) {
			LinkedList<T> pneq = pn.equivalents.get(pnOffset);
			if (pneq != null) {
				eq.addAll(pneq);
				for (T n : pneq) {
					weight += graph.getWeight(n);
				}
			}
		}
	}
	/**
	 * Iterates over each node in the given path
	 * @param it path
	 * @return sequence of primary nodes along the given path
	 */
	public static <T> Iterator<T> nodeIterator(Iterator<? extends PathNode<T>> it) {
		return nodeIterator(it, 0);
	}
	/**
	 * Iterates over each node in the given path
	 * @param it path
	 * @param offset number starting of nodes to skip 
	 * @return sequence of primary node along the given path
	 */
	public static <T> Iterator<T> nodeIterator(Iterator<? extends PathNode<T>> it, int offset) {
		Iterator<T> result = Iterators.concat(Iterators.transform(it, new Function<PathNode<T>, Iterator<T>>() {
			@Override
			public Iterator<T> apply(PathNode<T> input) {
				return input.path.iterator();
			}
		}));
		if (offset > 0) {
			Iterators.advance(result, offset);
		}
		return result;
	}
	/**
	 * Gets the node length of the given path
	 * @param it node path
	 * @return total node length
	 */
	public static <T> int nodeLength(Iterable<? extends PathNode<T>> it) {
		return nodeLength(it.iterator());
	}
	/**
	 * Gets the node length of the given path
	 * @param it node path
	 * @return total node length
	 */
	public static <T> int nodeLength(Iterator<? extends PathNode<T>> it) {
		int len = 0;
		while (it.hasNext()) {
			PathNode<T> pn = it.next();
			len += pn.length();
		}
		return len;
	}
	/**
	 * Gets the total weight of all nodes on the given path
	 * @param path path
	 * @return total weight
	 */
	public static <T> int totalWeight(Iterable<? extends PathNode<T>> it) {
		int weight = 0;
		for (PathNode<T> pn : it) {
			weight += pn.weight();
		}
		return weight;
	}
	/**
	 * Gets the total weight of all nodes on the given subpath
	 * @param it path
	 * @param offset number of starting nodes to skip
	 * @param length number of nodes to calculate weight of
	 * @return total weight
	 */
	public static <T> int totalWeight(Iterable<? extends PathNode<T>> it, int offset, int length, WeightedDirectedGraph<T> graph) {
		int weight = 0;
		for (PathNode<T> pn : it) {
			if (offset == 0 && pn.length() <= length) {
				// take entire node
				weight += pn.weight();
				length -= pn.length();
			} else {
				// pro-rata weight according to all the nodes in the positions of interest
				int toSkip = Math.min(pn.length(), offset);
				int toTake = Math.min(pn.length() - toSkip, length);
				for (int i = 0; i < toTake; i++) {
					weight += pn.weight(i + toSkip, graph);
				}
				offset -= toSkip;				
				length -= toTake;
			}
		}
		if (length != 0) throw new IllegalArgumentException("Too many nodes requested");
		if (offset != 0) throw new IllegalArgumentException("Skipping more nodes than exist on path");
		return weight;
	}
	public int indexOf(T node) {
		for (int i = 0; i < length(); i++) {
			if (path.get(i).equals(node)) {
				return i;
			}
		}
		if (equivalents != null) {
			for (int i = 0; i < length(); i++) {
				if (equivalents.get(i) != null && equivalents.get(i).contains(node)) {
					return i;
				} 
			}
		}
		return -1;
	}
	public boolean contains(T node) {
		return indexOf(node) >= 0;
	}
	@Override
	public String toString() {
		return toString(null);
	}
	public String toString(WeightedDirectedGraph<T> graph) {
		if (graph != null) {
			return String.format("[%4d]%s\t%s", nodeId, printAttributes(), graph.toString(path));
		} else {
			return String.format("[%4d]%s\t?", nodeId, printAttributes());
		}
	}
	protected String printAttributes() {
		return String.format(" l=%d\tn=%d\tw=%d", length(), nodes(), weight());
	}
}