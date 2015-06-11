package au.edu.wehi.idsv.graph;

import java.util.Iterator;

public class WeightedSequenceGraphNodeUtil {
	/**
	 * Gets the node length of the given path
	 * @param it node path
	 * @return total node length
	 */
	public static int nodeLength(Iterable<? extends WeightedSequenceGraphNode> it) {
		return nodeLength(it.iterator());
	}
	/**
	 * Gets the node length of the given path
	 * @param it node path
	 * @return total node length
	 */
	public static int nodeLength(Iterator<? extends WeightedSequenceGraphNode> it) {
		int len = 0;
		while (it.hasNext()) {
			WeightedSequenceGraphNode pn = it.next();
			len += pn.length();
		}
		return len;
	}
	/**
	 * Gets the total weight of all nodes on the given path
	 * @param path path
	 * @return total weight
	 */
	public static int totalWeight(Iterable<? extends WeightedSequenceGraphNode> it) {
		int weight = 0;
		for (WeightedSequenceGraphNode pn : it) {
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
	public static int totalWeight(Iterable<? extends WeightedSequenceGraphNode> it, int offset, int length) {
		int weight = 0;
		for (WeightedSequenceGraphNode pn : it) {
			if (offset == 0 && pn.length() <= length) {
				// take entire node
				weight += pn.weight();
				length -= pn.length();
			} else {
				// pro-rata weight according to all the nodes in the positions of interest
				int toSkip = Math.min(pn.length(), offset);
				int toTake = Math.min(pn.length() - toSkip, length);
				for (int i = 0; i < toTake; i++) {
					weight += pn.weight(i + toSkip);
				}
				offset -= toSkip;				
				length -= toTake;
			}
		}
		if (length != 0) throw new IllegalArgumentException("Too many nodes requested");
		if (offset != 0) throw new IllegalArgumentException("Skipping more nodes than exist on path");
		return weight;
	}
}
