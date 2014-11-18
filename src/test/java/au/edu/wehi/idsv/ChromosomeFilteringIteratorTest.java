package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import htsjdk.samtools.util.CloseableIterator;

import java.util.Iterator;
import java.util.List;

import org.junit.Test;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;


public class ChromosomeFilteringIteratorTest extends TestHelper {
	private static class TestCloseableIterator<T> implements CloseableIterator<T> {
		private final Iterator<T> it;
		public boolean isClosed = false;
		public TestCloseableIterator(final Iterator<T> it) {
			this.it = it;
		}
		@Override
		public void close() {
			isClosed = true;
		}
		@Override
		public boolean hasNext() {
			return it.hasNext();
		}
		@Override
		public T next() {
			return it.next();
		}
		@Override
		public void remove() {
			it.remove();
		}
	}
	@Test
	public void should_return_correct_chr_records() {
		TestCloseableIterator<DirectedEvidence> ci = new TestCloseableIterator<DirectedEvidence>(ImmutableList.<DirectedEvidence>of(
				E(0, 1, FWD),
				E(1, 1, FWD),
				E(1, 2, FWD),
				E(1, 3, FWD),
				E(2, 1, FWD)
			).iterator());
		ChromosomeFilteringIterator<DirectedEvidence> it = new ChromosomeFilteringIterator<DirectedEvidence>(ci, 1, true);
		List<DirectedEvidence> result = Lists.newArrayList(it);
		assertEquals(3, result.size());
		assertEquals(1, result.get(0).getBreakendSummary().start);
		assertEquals(2, result.get(1).getBreakendSummary().start);
		assertEquals(3, result.get(2).getBreakendSummary().start);
		assertTrue(ci.isClosed);
		it.close();
	}
	@Test
	public void should_close_as_soon_as_possible() {
		TestCloseableIterator<DirectedEvidence> ci = new TestCloseableIterator<DirectedEvidence>(ImmutableList.<DirectedEvidence>of(
				E(0, 1, FWD),
				E(1, 1, FWD),
				E(2, 1, FWD),
				E(2, 1, FWD),
				E(2, 1, FWD),
				E(2, 1, FWD)
			).iterator());
		ChromosomeFilteringIterator<DirectedEvidence> it = new ChromosomeFilteringIterator<DirectedEvidence>(ci, 1, true);
		it.next();
		assertFalse(ci.isClosed);
		it.hasNext();
		assertTrue(ci.isClosed);
		it.close();
	}
	@Test
	public void should_not_early_close_unsorted() {
		TestCloseableIterator<DirectedEvidence> ci = new TestCloseableIterator<DirectedEvidence>(ImmutableList.<DirectedEvidence>of(
				E(0, 1, FWD),
				E(1, 1, FWD),
				E(2, 1, FWD),
				E(2, 1, FWD),
				E(2, 1, FWD),
				E(1, 2, FWD),
				E(2, 1, FWD)
			).iterator());
		ChromosomeFilteringIterator<DirectedEvidence> it = new ChromosomeFilteringIterator<DirectedEvidence>(ci, 1, false);
		it.next();
		assertFalse(ci.isClosed);
		it.next();
		assertFalse(ci.isClosed);
		it.close();
	}
	@Test
	public void should_short_circuit_sorted() {
		ImmutableList<DirectedEvidence> x = ImmutableList.<DirectedEvidence>of(
				E(0, 1, FWD),
				E(1, 1, FWD),
				E(2, 1, FWD),
				E(2, 1, FWD),
				E(2, 1, FWD),
				E(1, 2, FWD),
				E(2, 1, FWD));
		assertEquals(1, Lists.newArrayList(new ChromosomeFilteringIterator<DirectedEvidence>(x.iterator(), 1, true)).size());
		assertEquals(2, Lists.newArrayList(new ChromosomeFilteringIterator<DirectedEvidence>(x.iterator(), 1, false)).size());
	}
}
