package au.edu.wehi.idsv.util;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.util.CloseableIterator;

import java.util.List;

import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Ordering;


public class AutoClosingMergedIteratorTest extends TestHelper {
	public static class CountingIterator extends AbstractIterator<Integer> implements CloseableIterator<Integer> {
		private int count;
		private int current = 0; 
		public volatile boolean isClosed = false;
		public CountingIterator(int count) { this.count = count; }
		@Override
		public void close() {
			isClosed = true;
		}
		@Override
		protected Integer computeNext() {
			if (++current <= count) return current;
			return endOfData();
		}
	}
	@Test
	public void should_terminate_async_reader_thread_for_underlying_stream() throws InterruptedException {
		List<SAMRecord> list = ImmutableList.of(Read(0, 1, 1), Read(0, 1, 1), Read(0, 1, 1), Read(0, 1, 1), Read(0, 2, 1));
		AsyncBufferedIterator<SAMRecord> abi = new AsyncBufferedIterator<SAMRecord>(list.iterator(), 1, 1);
		AutoClosingMergedIterator<SAMRecord> merged = new AutoClosingMergedIterator<SAMRecord>(ImmutableList.of(
				abi), new SAMRecordCoordinateComparator());
		merged.next();
		Thread.sleep(50);
		assertNotNull(AsyncBufferedIteratorTest.getThreadWithName(abi.getBackgroundThreadName()));
		merged.close();
		Thread.sleep(50);
		assertNull(AsyncBufferedIteratorTest.getThreadWithName(abi.getBackgroundThreadName()));
	}
	@Test
	public void close_should_close_all() {
		CountingIterator it1 = new CountingIterator(1);
		CountingIterator it2 = new CountingIterator(8);
		CountingIterator it3 = new CountingIterator(3);
		CountingIterator it4 = new CountingIterator(10);
		CountingIterator it5 = new CountingIterator(15);
		AutoClosingMergedIterator<Integer> merged = new AutoClosingMergedIterator<Integer>(ImmutableList.of(
				it1, it2, it3, it4, it5), Ordering.natural());
		// pull 10 elements
		for (int i = 0; i < 5; i++) assertEquals(1, (int)merged.next());
		assertTrue(it1.isClosed); // run out of elements
		assertFalse(it2.isClosed);
		for (int i = 0; i < 4; i++) assertEquals(2, (int)merged.next());
		for (int i = 0; i < 4; i++) assertEquals(3, (int)merged.next());
		assertTrue(it3.isClosed);
		assertFalse(it2.isClosed);
		merged.close();
		assertTrue(it1.isClosed);
		assertTrue(it2.isClosed);
		assertTrue(it3.isClosed);
		assertTrue(it4.isClosed);
		assertTrue(it5.isClosed);
	}
}
