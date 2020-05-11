package au.edu.wehi.idsv.util;

import com.google.common.collect.ImmutableList;

import java.util.concurrent.*;


public class ThreadPoolTest {
	public class TestCallable implements Callable<Void> {
		@Override
		public Void call() throws InterruptedException {
			Thread.sleep(500);
			return null;
		}
	}
	//@Test
	public void callable_void_should_return() throws InterruptedException {
		ExecutorService	es = Executors.newFixedThreadPool(2);
		es.invokeAll(ImmutableList.of(new TestCallable(), new TestCallable(), new TestCallable(), new TestCallable(), new TestCallable(), new TestCallable()));
	}
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
	//@Test // locks up since the blocking thread is in the thread pool using up the only thread
	public void blocking_tasks_cannot_be_submitted_from_a_worker_thread() throws InterruptedException {
		ExecutorService	es = Executors.newSingleThreadExecutor();
		es.invokeAll(ImmutableList.of(new RecursiveTask(es, 0)));
	}
	public class RecursiveForkJoinAction extends RecursiveAction {
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
	//@Test// yet breakend assembly deadlocks (maybe a resource threading issue?)
	public void blocking_tasks_can_be_submitted_from_a_fork_join_task() throws InterruptedException {
		ForkJoinPool p = new ForkJoinPool(1);
		p.invoke(new RecursiveForkJoinAction(0));
	}
}
