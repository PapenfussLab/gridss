package au.edu.wehi.idsv.util;

import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveAction;

import org.junit.Ignore;
import org.junit.Test;

import com.google.common.collect.ImmutableList;

public class ThreadPoolTest {
	public class RecursiveTask implements Callable<Integer> {
		private ExecutorService threadPool;
		private int depth;
		public RecursiveTask(ExecutorService threadPool, int depth) {
			this.threadPool = threadPool;
			this.depth = depth;
		}
		@Override
		public Integer call() throws Exception {
			if (depth < 8) {
				threadPool.invokeAll(ImmutableList.of(new RecursiveTask(threadPool, depth + 1)));
			}
			return depth;
		}
	}
	//@Test
	@Ignore // locks up since the blocking thread is in the thread pool using up the only thread
	public void blocking_tasks_cannot_be_submitted_from_a_worker_thread() throws InterruptedException {
		ExecutorService	es = Executors.newSingleThreadExecutor();
		es.invokeAll(ImmutableList.of(new RecursiveTask(es, 0)));
	}
	public class RecursiveForkJoinAction extends RecursiveAction {
		/** */
		private static final long serialVersionUID = 1299799042074340392L;
		private int depth;
		public RecursiveForkJoinAction(int depth) {
			this.depth = depth;
		}
		@Override
		protected void compute() {
			if (depth < 8) {
				super.
				invokeAll(new RecursiveForkJoinAction(depth + 1), new RecursiveForkJoinAction(depth + 1));
			}
		}
	}
	@Test
	public void blocking_tasks_can_be_submitted_from_a_fork_join_task() throws InterruptedException {
		ForkJoinPool p = new ForkJoinPool(1);
		p.invoke(new RecursiveForkJoinAction(0));
	}
}
