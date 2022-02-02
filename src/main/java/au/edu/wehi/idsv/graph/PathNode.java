package au.edu.wehi.idsv.graph;

import com.google.common.base.Function;
import com.google.common.collect.Iterators;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

public class PathNode<T> implements WeightedSequenceGraphNode {
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
	protected final WeightedDirectedGraph<T> graph;
	/**
	 * Creates a new path graph node from the given node sequence
	 * @param nodes node sequence
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
		this.graph = graph;
	}
	/**
	 * Merges the given path into the current path node
	 * @param path node sequence
	 * @param graph parent graph
	 * @return 
	 */
	public static <T> PathNode<T> concat(Iterable<? extends PathNode<T>> path, WeightedDirectedGraph<T> graph) {
		int length = WeightedSequenceGraphNodeUtil.nodeLength(path);
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
		this.path = new ArrayList<T>(length);
		int offset = 0;
		for (PathNode<T> pn : nodes) {
			for (int j = 0; j < pn.length() && offset - startOffset < length; j++, offset++) {
				if (offset < startOffset) {
					// skip
				} else {
					T n = pn.path.get(j);
					path.add(n);
				}
			}
		}
		this.graph = graph;
		recalculateWeight();
	}
	private void recalculateWeight() {
		weight = 0;
		for (T n : path) {
			weight += graph.getWeight(n);
		}
	}
	public List<T> getPath() { return path; }
	/**
	 * Gets node paths that have been merged into this path
	 * @return
	 */
	public T first() { return path.get(0); }
	public T last() { return path.get(length() - 1); }
	public int length() { return path.size(); }
	public int nodes() {
		int nodeCount = path.size();
		return nodeCount;
	}
	public int weight() {
		return weight;
	}
	public int weight(int offset) {
		int w = graph.getWeight(path.get(offset));
		return w;
	}
	public int maxNodeWeight() {
		int maxNodeWeight = 0;
		for (T n : getPath()) {
			maxNodeWeight = Math.max(maxNodeWeight, graph.getWeight(n));
		}
		return maxNodeWeight;
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
	public int indexOf(T node) {
		for (int i = 0; i < length(); i++) {
			if (path.get(i).equals(node)) {
				return i;
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