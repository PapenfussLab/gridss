package au.edu.wehi.idsv.graph;

import java.util.ArrayDeque;
import java.util.Iterator;
import java.util.Queue;

import com.google.common.collect.AbstractIterator;

/**
 * Streaming maximal clique caller wrapper
 * 
 * @author Daniel Cameron
 */
public class RectangleGraphMaximalCliqueIterator extends AbstractIterator<GraphNode> {
	private final Queue<GraphNode> buffer = new ArrayDeque<GraphNode>();
	private RectangleGraphMaximalCliqueCalculator calc = new RectangleGraphMaximalCliqueCalculator();
	private Iterator<GraphNode> it;
	public RectangleGraphMaximalCliqueIterator(Iterator<GraphNode> it) {
		this.it = it;
	}
	@Override
	protected GraphNode computeNext() {
		while (buffer.isEmpty() && it.hasNext()) {
			buffer.addAll(calc.next(it.next()));
		}
		if (buffer.isEmpty() && calc != null) {
			buffer.addAll(calc.complete());
			calc = null;
		}
		while (!buffer.isEmpty()) {
			return buffer.poll();
		}
		return endOfData();
	}
}