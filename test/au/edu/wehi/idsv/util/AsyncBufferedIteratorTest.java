package au.edu.wehi.idsv.util;

import static org.junit.Assert.*;

import java.util.List;

import htsjdk.samtools.util.CloseableIterator;

import org.junit.Test;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import com.google.common.primitives.Ints;


public class AsyncBufferedIteratorTest {
	public class CIT extends AbstractIterator<Integer> implements CloseableIterator<Integer> {
		private int count;
		public volatile boolean isClosed = false;
		public CIT(int count) { this.count = count; }
		@Override
		public void close() {
			isClosed = true;
		}
		@Override
		protected Integer computeNext() {
			if (count > 0) return count--;
			return endOfData();
		}
	}
	@Test
	public void should_return_underlying_records() {
		List<Integer> list = Ints.asList(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12);
		assertArrayEquals(Ints.toArray(list), Ints.toArray(Lists.newArrayList(new AsyncBufferedIterator<Integer>(list.iterator(), 1, 1))));
		assertEquals(1024, Iterators.size(new AsyncBufferedIterator<Integer>(new CIT(1024), 8, 8)));
	}
	@Test
	public void should_close_underlying() {
		CIT it = new CIT(5);
		AsyncBufferedIterator<Integer> abi = new AsyncBufferedIterator<Integer>(it, 1, 1);
		abi.close();
		assertTrue(it.isClosed);
	}
	@Test
	public void should_close_when_underlying_end_of_stream_reached() throws InterruptedException {
		// 0 1 2 3
		// 4 5 6 7
		// 8 91011
		// 1213
		CIT it = new CIT(13);
		AsyncBufferedIterator<Integer> abi = new AsyncBufferedIterator<Integer>(it, 4, 2);
		abi.next(); assertFalse(it.isClosed); // first 4 records removed from buffer 
		abi.next(); assertFalse(it.isClosed);
		abi.next(); assertFalse(it.isClosed);
		abi.next(); assertFalse(it.isClosed);
		abi.next(); // next 4 records removed
		Thread.sleep(1); // let the consumer thread fill up the newly opened buffer slot with the final two records
		assertTrue(it.isClosed); // should have now closed
		abi.close();
	}
}
