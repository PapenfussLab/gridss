package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import htsjdk.samtools.util.CloseableIterator;

import java.io.File;
import java.util.Iterator;
import java.util.List;

import org.junit.Test;

import au.edu.wehi.idsv.metrics.RelevantMetrics;

import com.google.common.collect.Lists;


public class EvidenceSourceTest extends IntermediateFilesTest {
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
	public class TestEvidenceSource extends EvidenceSource {
		public TestEvidenceSource(boolean perChr, File input) {
			super(getCommandlineContext(perChr), input);
		}
		private boolean processed = false;
		@Override
		protected boolean isProcessingComplete() { return processed; }
		@Override
		protected void process() { processed = true; }
		@Override
		protected CloseableIterator<DirectedEvidence> perChrIterator(String chr) {
			int index = processContext.getDictionary().getSequence(chr).getSequenceIndex();
			List<DirectedEvidence> list = Lists.newArrayList();
			list.add(SCE(BWD, Read(index, 1, "50S50M")));
			list.add(SCE(BWD, Read(index, 2, "50S50M")));
			list.add(SCE(BWD, Read(index, 3, "50S50M")));
			lastIt = new TestCloseableIterator<DirectedEvidence>(list.iterator());
			return lastIt;
		}
		@Override
		protected CloseableIterator<DirectedEvidence> singleFileIterator() {
			List<DirectedEvidence> list = Lists.newArrayList();
			list.add(SCE(BWD, Read(0, 1, "50S50M")));
			list.add(SCE(BWD, Read(0, 2, "50S50M")));
			list.add(SCE(BWD, Read(0, 3, "50S50M")));
			list.add(SCE(BWD, Read(1, 4, "50S50M")));
			list.add(SCE(BWD, Read(1, 5, "50S50M")));
			list.add(SCE(BWD, Read(1, 6, "50S50M")));
			list.add(SCE(BWD, Read(2, 7, "50S50M")));
			list.add(SCE(BWD, Read(2, 8, "50S50M")));
			list.add(SCE(BWD, Read(2, 9, "50S50M")));
			lastIt = new TestCloseableIterator<DirectedEvidence>(list.iterator());
			return lastIt;
		}
		public TestCloseableIterator<DirectedEvidence> lastIt = null;
		@Override
		public RelevantMetrics getMetrics() {
			return null;
		}
	}
	@Test
	public void getFileIntermediateDirectoryBasedOn_should_be_input_file() {
		EvidenceSource es = new TestEvidenceSource(false, input);
		assertEquals(input, es.getFileIntermediateDirectoryBasedOn());
	}
	@Test
	public void CloseableChromosomeFilterIterator_should_return_correct_chr_records() {
		TestEvidenceSource es = new TestEvidenceSource(false, input);
		List<DirectedEvidence> result = Lists.newArrayList(es.iterator("polyACGT"));
		assertEquals(3, result.size());
		assertEquals(4, result.get(0).getBreakendSummary().start);
		assertEquals(5, result.get(1).getBreakendSummary().start);
		assertEquals(6, result.get(2).getBreakendSummary().start);
	}
	@Test
	public void CloseableChromosomeFilterIterator_should_close_as_soon_as_possible() {
		TestEvidenceSource es = new TestEvidenceSource(false, input);
		Iterator<DirectedEvidence> it = es.iterator("polyACGT");
		assertEquals(1, it.next().getBreakendSummary().referenceIndex);
		assertEquals(1, it.next().getBreakendSummary().referenceIndex);
		assertEquals(1, it.next().getBreakendSummary().referenceIndex);
		assertFalse(es.lastIt.isClosed);
		it.hasNext();
		assertTrue(es.lastIt.isClosed);
	}
	@Test
	public void PerChromosomeAggregateIterator_should_return_all_chr() {
		TestEvidenceSource es = new TestEvidenceSource(true, input);
		List<DirectedEvidence> result = Lists.newArrayList(es.iterator());
		assertEquals(3 * getSequenceDictionary().size(), result.size());
	}
	@Test
	public void PerChromosomeAggregateIterator_should_close_as_soon_as_possible() {
		TestEvidenceSource es = new TestEvidenceSource(true, input);
		Iterator<DirectedEvidence> it = es.iterator();
		it.next();
		it.next();
		it.next();
		TestCloseableIterator<DirectedEvidence> underlyingIt = es.lastIt;
		assertFalse(underlyingIt.isClosed);
		it.hasNext(); // we hit the end of our current iterator
		assertTrue(underlyingIt.isClosed);
		
		es = new TestEvidenceSource(true, input);
		it = es.iterator();
		it.next();
		it.next();
		it.next();
		underlyingIt = es.lastIt;
		assertFalse(underlyingIt.isClosed);
		it.next(); // next() should also close
		assertTrue(underlyingIt.isClosed);
	}
}
