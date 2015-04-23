package au.edu.wehi.idsv.graph;

import java.util.Comparator;
import java.util.Iterator;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

/**
 * Merges GraphNodes with the same coordinates
 * @author Daniel Cameron
 *
 */
public class RectangleGraphNodeMergingIterator extends AbstractIterator<RectangleGraphNode> {
	private final PeekingIterator<RectangleGraphNode> it;
	private final Comparator<RectangleGraphNode> expectedOrder;
	public RectangleGraphNodeMergingIterator(Comparator<RectangleGraphNode> expectedOrder, Iterator<RectangleGraphNode> it) {
		this.it = Iterators.peekingIterator(it);
		this.expectedOrder = expectedOrder;
	}
	@Override
	protected RectangleGraphNode computeNext() {
		if (it.hasNext()) {
			RectangleGraphNode node = it.next();
			while (it.hasNext() && it.peek().isSameCoordinate(node)) {
				node = new RectangleGraphNode(node.startX, node.endX, node.startY, node.endY, node.weight + it.next().weight);
			}
			if (expectedOrder != null) {
				if (it.hasNext() && expectedOrder.compare(node, it.peek()) > 0) {
					throw new IllegalStateException(String.format("Unexpected out of order sequence: %s before %s", node, it.peek()));
				}
			}
			return node;
		}
		return endOfData();
	}
}
