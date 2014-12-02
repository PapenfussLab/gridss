package au.edu.wehi.idsv.graph;

import java.util.Iterator;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

/**
 * Merges GraphNodes with the same coordinates
 * @author Daniel Cameron
 *
 */
public class GraphNodeMergingIterator extends AbstractIterator<GraphNode> {
	private final PeekingIterator<GraphNode> it;
	public GraphNodeMergingIterator(Iterator<GraphNode> it) {
		this.it = Iterators.peekingIterator(it);
	}
	@Override
	protected GraphNode computeNext() {
		if (it.hasNext()) {
			GraphNode node = it.next();
			while (it.hasNext() && it.peek().isSameCoordinate(node)) {
				node = new GraphNode(node.startX, node.endX, node.startY, node.endY, node.weight + it.next().weight);
			}
			return node;
		}
		return endOfData();
	}
}
