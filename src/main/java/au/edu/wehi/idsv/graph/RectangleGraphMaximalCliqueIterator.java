package au.edu.wehi.idsv.graph;

import java.util.ArrayDeque;
import java.util.Collection;
import java.util.Iterator;
import java.util.Queue;

import au.edu.wehi.idsv.visualisation.TrackedState;
import com.google.common.collect.AbstractIterator;

/**
 * Streaming maximal clique caller wrapper
 * 
 * @author Daniel Cameron
 */
public class RectangleGraphMaximalCliqueIterator extends AbstractIterator<RectangleGraphNode> implements TrackedState {
	private final Queue<RectangleGraphNode> buffer = new ArrayDeque<RectangleGraphNode>();
	private RectangleGraphMaximalCliqueCalculator calc = new RectangleGraphMaximalCliqueCalculator();
	private Iterator<RectangleGraphNode> it;
	public RectangleGraphMaximalCliqueIterator(Iterator<RectangleGraphNode> it) {
		this.it = it;
	}
	@Override
	protected RectangleGraphNode computeNext() {
		while (buffer.isEmpty() && it.hasNext()) {
			RectangleGraphNode nextEvidence = it.next(); 
			buffer.addAll(calc.next(nextEvidence));
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

	@Override
	public String[] trackedNames() {
		return calc.trackedNames();
	}

	@Override
	public Object[] trackedState() {
		return calc.trackedState();
	}

	@Override
	public Collection<TrackedState> trackedObjects() {
		calc.trackedObjects();
	}
}