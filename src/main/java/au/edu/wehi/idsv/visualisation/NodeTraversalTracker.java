package au.edu.wehi.idsv.visualisation;

public interface NodeTraversalTracker<T> {
	/**
	 * Track the visitation of the given node
	 * @param node
	 */
	void track(T node);
	/**
	 * Count of nodes tracked
	 * @return
	 */
	int count();
}
