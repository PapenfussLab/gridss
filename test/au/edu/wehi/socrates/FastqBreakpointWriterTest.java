package au.edu.wehi.socrates;

import static org.junit.Assert.*;

import java.util.List;

import org.junit.Test;

import com.google.common.collect.Lists;

import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;

public class FastqBreakpointWriterTest extends TestHelper {
	public class MockWriter implements FastqWriter {
		public boolean closeCalled = false;
		public List<FastqRecord> list = Lists.newArrayList();
		@Override
		public void write(FastqRecord rec) {
			list.add(rec);
		}
		@Override
		public void close() {
			closeCalled = true;
		}
	}
	@Test
	public void should_call_close() {
		MockWriter w = new MockWriter();
		FastqBreakpointWriter fw = new FastqBreakpointWriter(w);
		fw.close();
		assertTrue(w.closeCalled);
	}
	public void should_call_rec() {
		MockWriter w = new MockWriter();
		FastqBreakpointWriter fw = new FastqBreakpointWriter(w);
		fw.write(new SoftClipEvidence(getContext(), BreakpointDirection.Forward, Read(0, 1, "1M1S")));
		assertEquals(1, w.list.size());
	}
}
