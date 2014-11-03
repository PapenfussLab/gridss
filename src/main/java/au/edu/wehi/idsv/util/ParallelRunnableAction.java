package au.edu.wehi.idsv.util;

import java.util.Collection;
import java.util.List;
import java.util.concurrent.ForkJoinTask;
import java.util.concurrent.RecursiveAction;

import com.google.common.base.Function;
import com.google.common.collect.Lists;

/**
 * Wrapper utility class as ForkJoinPool does not allow invokeAll on a ForkJoinTask collection
 * @author cameron.d
 * @param <T>
 *
 */
public class ParallelRunnableAction extends RecursiveAction {
	/**
	 * 
	 */
	private static final long serialVersionUID = 299198590877452953L;
	private List<Runnable> tasks;
	public ParallelRunnableAction(Collection<? extends Runnable> tasks) {
		this.tasks = Lists.newArrayList(tasks);
	}
	@SuppressWarnings("rawtypes")
	@Override
	protected void compute() {
		invokeAll(Lists.<Runnable, ForkJoinTask>transform(tasks, new Function<Runnable, ForkJoinTask>() {
			@Override
			public ForkJoinTask apply(Runnable input) {
				return ForkJoinTask.adapt(input);
			}
		}));
	}
}
