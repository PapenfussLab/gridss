package au.edu.wehi.idsv.util;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Iterator;
import java.util.List;

import org.junit.Test;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;


public class DuplicatingIterableTest {
	@Test
	public void should_return_underlying_iterator() {
		List<Integer> list = ImmutableList.of(0, 1, 2, 3, 4, 5, 6, 7);
		assertEquals(list, Lists.newArrayList(new DuplicatingIterable<Integer>(list.iterator(), 1).iterator()));
	}
	@Test
	public void should_allow_traversal_difference_up_to_and_including_buffer_size() {
		List<Integer> list = ImmutableList.of(0, 1, 2, 3, 4, 5, 6, 7);
		DuplicatingIterable<Integer> dib = new DuplicatingIterable<Integer>(list.iterator(), 2);
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
	public void should_clone_from_underlying_position() {
		List<Integer> list = ImmutableList.of(0, 1, 2, 3, 4, 5, 6, 7);
		DuplicatingIterable<Integer> dib = new DuplicatingIterable<Integer>(list.iterator(), 2);
		Iterator<Integer> it1 = dib.iterator();
		it1.next();
		assertEquals(1, (int)dib.iterator().next());
	}
	@Test
	public void should_block_when_reading_too_far_ahead_of_other_iterators() throws InterruptedException {
		List<Integer> list = ImmutableList.of(0, 1, 2, 3, 4, 5, 6, 7);
		DuplicatingIterable<Integer> dib = new DuplicatingIterable<Integer>(list.iterator(), 2);
		Iterator<Integer> it1 = dib.iterator();
		ConsumerThread<Integer> t = new ConsumerThread<Integer>(dib.iterator(), 3);
		t.start();
		Thread.sleep(1);
		assertTrue(t.isAlive());
		it1.next(); // consume record thus allowing the thread to read the third record
		Thread.sleep(1);
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
}