package au.edu.wehi.idsv.util;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.collect.PeekingIterator;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;

import static org.junit.Assert.*;


public class DuplicatingIterableTest {
	@Test
	public void should_return_underlying_iterator() {
		List<Integer> list = ImmutableList.of(0, 1, 2, 3, 4, 5, 6, 7);
		assertEquals(list, Lists.newArrayList(new DuplicatingIterable<Integer>(1, list.iterator(), 1).iterator()));
	}
	@Test
	public void should_allow_nonblocking_traversal_difference_up_to_and_including_buffer_size() {
		List<Integer> list = ImmutableList.of(0, 1, 2, 3, 4, 5, 6, 7);
		DuplicatingIterable<Integer> dib = new DuplicatingIterable<Integer>(2, list.iterator(), 2);
		Iterator<Integer> it1 = dib.iterator();
		Iterator<Integer> it2 = dib.iterator();
		it1.next();
		it1.next();
		// it1 is now 2 ahead of it2
		it2.next();
		it2.next();
		it2.next();
		it2.next();
		// it2 is now 2 ahead of it1
	}
	@Test
	public void should_iterate_from_start() {
		List<Integer> list = ImmutableList.of(0, 1, 2, 3, 4, 5, 6, 7);
		DuplicatingIterable<Integer> dib = new DuplicatingIterable<Integer>(2, list.iterator(), 2);
		Iterator<Integer> it1 = dib.iterator();
		it1.next();
		assertEquals(0, (int)dib.iterator().next());
	}
	@Test
	public void should_block_when_reading_too_far_ahead_of_other_iterators() throws InterruptedException {
		List<Integer> list = ImmutableList.of(0, 1, 2, 3, 4, 5, 6, 7);
		DuplicatingIterable<Integer> dib = new DuplicatingIterable<Integer>(2, list.iterator(), 2);
		Iterator<Integer> it1 = dib.iterator();
		ConsumerThread<Integer> t = new ConsumerThread<Integer>(dib.iterator(), 3);
		t.start();
		Thread.sleep(10);
		assertTrue(t.isAlive());
		it1.next(); // consume record thus allowing the thread to read the third record
		for (int i = 0; i < 128; i++) {
			if (t.isAlive()) {
				Thread.sleep(1);
			}
		}
		assertFalse(t.isAlive());
	}
	public class ConsumerThread<T> extends Thread {
		private Iterator<T> it;
		private int recordsToConsume;
		public ConsumerThread(Iterator<T> it, int recordsToConsume) {
			this.it = it;
			this.recordsToConsume = recordsToConsume;
		}
		@Override
		public void run() {
			for (int i = 0; i < recordsToConsume; i++) {
				it.next();
			}
		}
	}
	@Test
	public void should_allow_parallel_duplicating_producer_pattern() throws InterruptedException {
		List<Integer> list = new ArrayList<Integer>();
		for (int i = 0; i < 1024; i++) {
			list.add(i);
		}
		destination = new ArrayBlockingQueue<Integer>(1);
		int threads = 8;
		DuplicatingIterable<Integer> dib = new DuplicatingIterable<Integer>(threads, list.iterator(), 1);
		Thread[] producers = new Thread[threads];
		for (int i = 0; i < producers.length; i++) {
			producers[i] = new ProducerThread(dib.iterator());
		}
		for (int i = 0; i < producers.length; i++) {
			producers[i].start();
		}
		for (int i = 0; i < list.size() * producers.length; i++) {
			Integer result = destination.take();
			assert(result <= i); // no guarantees about threading order
		}
		Thread.sleep(50);
		assertTrue(destination.isEmpty());
		Thread.sleep(50);
		for (int i = 0; i < producers.length; i++) {
			assertFalse(producers[i].isAlive());
		}
	}
	private ArrayBlockingQueue<Integer> destination;
	public class ProducerThread extends Thread {
		private Iterator<Integer> it;
		public ProducerThread(Iterator<Integer> it) {
			this.it = it;
		}
		@Override
		public void run() {
			while (it.hasNext()) {
				try {
					destination.put(it.next());
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
	}
	@Test
	public void should_raise_on_all_iterators() {
		int exceptionsFound = 0;
		int n = 4;
		DuplicatingIterable<Integer> di = new DuplicatingIterable<Integer>(n, new ErrorIterator<>(), 16);
		for (int i = 0; i < n; i++) {
			PeekingIterator<Integer> it = di.iterator();
			try {
				it.next();
			} catch (RuntimeException e) {
				exceptionsFound++;
			}
		}
		assertEquals(n, exceptionsFound);
	}
}