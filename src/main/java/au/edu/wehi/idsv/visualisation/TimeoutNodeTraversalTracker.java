package au.edu.wehi.idsv.visualisation;

import au.edu.wehi.idsv.util.AlgorithmRuntimeSafetyLimitExceededException;

public class TimeoutNodeTraversalTracker<T> implements NodeTraversalTracker<T> {
	public TimeoutNodeTraversalTracker(int maxTraversalCount) {
		this.maxTraversalCount = maxTraversalCount;
	}
	private int maxTraversalCount;
	private int nodeTaversalCount;
	@Override
	public void track(T node) {
		if (++nodeTaversalCount >= maxTraversalCount) {
			throw new AlgorithmRuntimeSafetyLimitExceededException(String.format(
				"Path collapse not complete after %d node traversals. Aborting path collapse whilst processing \"%s\".",
					nodeTaversalCount,
					node.toString()));
		}
	}
	public int count() {
		return nodeTaversalCount;
	}

}
