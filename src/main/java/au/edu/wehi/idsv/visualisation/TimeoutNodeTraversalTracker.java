package au.edu.wehi.idsv.visualisation;

import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.util.AlgorithmRuntimeSafetyLimitExceededException;

public class TimeoutNodeTraversalTracker<T> implements NodeTraversalTracker<T> {
	private int nodeTaversalCount;
	@Override
	public void track(T node) {
		if (++nodeTaversalCount >= Defaults.COLLAPSE_PATH_MAX_TRAVERSAL) {
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
