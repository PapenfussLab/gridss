package au.edu.wehi.idsv.util;

import java.util.Collection;
import java.util.concurrent.ForkJoinTask;
import java.util.concurrent.RecursiveAction;

/**
 * Wrapper utility class as ForkJoinPool does not allow invokeAll on a ForkJoinTask collection
 * @author cameron.d
 * @param <T>
 *
 */
public class ParallelAction<T extends ForkJoinTask<?>> extends RecursiveAction {
	/**
	 * 
	 */
	private static final long serialVersionUID = 299198590877452953L;
	private Collection<T> tasks;
	public ParallelAction(Collection<T> tasks) {
		this.tasks = tasks;
	}
	@Override
	protected void compute() {
		invokeAll(tasks);
	}

}
