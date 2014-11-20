package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.junit.Test;

import com.google.common.base.Function;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;


public class PerChromosomeAggregateIteratorTest extends TestHelper {
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
	public void PerChromosomeAggregateIterator_should_call_and_concatenate_iterator_generator_function_for_each_chr_in_dictionary_order() {
		SAMSequenceDictionary dict = getContext().getDictionary();
		PerChromosomeAggregateIterator<String> it = new PerChromosomeAggregateIterator<String>(dict, new Function<String, Iterator<String>>() {
			@Override
			public Iterator<String> apply(String input) {
				return ImmutableList.of(input).iterator();
			}
		});
		List<String> result = Lists.newArrayList(it);
		assertEquals(dict.getSequences().size(), result.size());
		for (int i = 0; i < result.size(); i++) {
			assertEquals(dict.getSequence(i).getSequenceName(), result.get(i));
		}
	}
	@Test
	public void PerChromosomeAggregateIterator_should_close_as_soon_as_possible() {
		SAMSequenceDictionary dict = getContext().getDictionary();
		final List<TestCloseableIterator<String>> createdIt = new ArrayList<>();
		PerChromosomeAggregateIterator<String> it = new PerChromosomeAggregateIterator<String>(dict, new Function<String, Iterator<String>>() {
			@Override
			public Iterator<String> apply(String input) {
				TestCloseableIterator<String> i = new TestCloseableIterator<String>(ImmutableList.of(input, input).iterator());
				createdIt.add(i);
				return i;
			}
		});
		assertEquals(0, createdIt.size()); // No iterator exists yets
		it.hasNext();
		assertFalse(createdIt.get(0).isClosed); // first iterator created
		it.next();
		assertFalse(createdIt.get(0).isClosed);
		it.next();
		// exhausted but we don't know it yet
		it.hasNext();
		assertTrue(createdIt.get(0).isClosed);
		it.close();
	}
}
