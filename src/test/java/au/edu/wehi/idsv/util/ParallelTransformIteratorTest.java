package au.edu.wehi.idsv.util;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.function.Function;

import org.junit.Test;

import au.edu.wehi.idsv.util.AsyncBufferedIteratorTest.CIT;

import com.google.common.collect.Lists;
import com.google.common.primitives.Ints;


public class ParallelTransformIteratorTest {
	@Test
	public void should_apply_transform() {
		for (int i = 1; i < 4; i++) {
			List<Integer> list = Ints.asList(0, 1, 2, 3);
			Function<Integer, Integer> f = n -> n + 1;
			ParallelTransformIterator<Integer, Integer> it = new ParallelTransformIterator<Integer, Integer>(list.iterator(), f, i, Runnable::run);
			List<Integer> results = Lists.newArrayList(it);
			assertEquals(results, Ints.asList(1, 2, 3, 4));
		}
	}
	@Test
	public void lookahead_should_determine_number_of_records_to_read_ahead() {
		CIT cit = new CIT(16);
		ParallelTransformIterator<Integer, Integer> it = new ParallelTransformIterator<Integer, Integer>(cit, n -> n, 4, Runnable::run);
		it.next();
		assertEquals(16-4-1, cit.recordsleft);
	}
	@Test
	public void should_not_start_iteration_until_next_is_called() {
		CIT cit = new CIT(16);
		ParallelTransformIterator<Integer, Integer> it = new ParallelTransformIterator<Integer, Integer>(cit, n -> n, 4, Runnable::run);
		assertEquals(16, cit.recordsleft);
		it.next();
		for (int i = 1; i < 16; i++) it.next();
		assertEquals(0, cit.recordsleft);
	}
	@Test
	public void should_retain_iteration_order() {
		ExecutorService threadpool = Executors.newFixedThreadPool(4);
		CIT cit = new CIT(32);
		ParallelTransformIterator<Integer, Integer> it = new ParallelTransformIterator<Integer, Integer>(cit, n -> {
			try {
				Thread.sleep(n);
			} catch (InterruptedException e) {
			}
			return n;
		}, 4, threadpool);
		for (int i = 32; i > 0; i--) assertEquals(i, (int)it.next());
	}
}
