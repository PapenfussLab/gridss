package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;

import java.util.List;

import org.junit.Test;

import com.google.common.collect.Lists;

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
	@Test
	public void should_call_rec() {
		MockWriter w = new MockWriter();
		FastqBreakpointWriter fw = new FastqBreakpointWriter(w);
		fw.write(new SoftClipEvidence(SES(), BreakendDirection.Forward, Read(0, 1, "1M1S")));
		assertEquals(1, w.list.size());
		fw.close();
	}
	@Test
	public void write_breakend_soft_clips() {
		MockWriter w = new MockWriter();
		FastqBreakpointWriter fw = new FastqBreakpointWriter(w);
		SoftClipEvidence e = new SoftClipEvidence(SES(), BreakendDirection.Forward, Read(0, 1, "10M10S"));
		fw.write(e);
		SAMRecord realigned = new SAMRecord(getContext().getBasicSamHeader());
		realigned.setReadName(w.list.get(0).getReadHeader());
		realigned.setReadBases(B("AAACCCGGGT"));
		realigned.setCigarString("1S5M4S");
		realigned.setBaseQualities(new byte[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9} );		
		fw.writeRealignment(realigned, 1);
		assertEquals(3, w.list.size());
		assertEquals("A", w.list.get(1).getReadString());
		assertEquals(0, BreakpointFastqEncoding.getEncodedBreakendOffset(w.list.get(1).getReadHeader()));
		assertEquals("GGGT", w.list.get(2).getReadString());
		assertEquals(6, BreakpointFastqEncoding.getEncodedBreakendOffset(w.list.get(2).getReadHeader()));
		fw.close();
	}
	@Test
	public void should_consider_soft_clip_length() {
		MockWriter w = new MockWriter();
		FastqBreakpointWriter fw = new FastqBreakpointWriter(w);
		SoftClipEvidence e = new SoftClipEvidence(SES(), BreakendDirection.Forward, Read(0, 1, "10M10S"));
		fw.write(e);
		SAMRecord realigned = new SAMRecord(getContext().getBasicSamHeader());
		realigned.setReadName(w.list.get(0).getReadHeader());
		realigned.setReadBases(B("AAACCCGGGT"));
		realigned.setCigarString("1S5M4S");
		realigned.setBaseQualities(new byte[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9} );		
		fw.writeRealignment(realigned, 2);
		assertEquals(2, w.list.size());
		fw.close();
	}
	@Test
	public void should_consider_realigned_orientiation() {
		MockWriter w = new MockWriter();
		FastqBreakpointWriter fw = new FastqBreakpointWriter(w);
		SoftClipEvidence e = new SoftClipEvidence(SES(), BreakendDirection.Forward, Read(0, 1, "10M10S"));
		fw.write(e);
		SAMRecord realigned = new SAMRecord(getContext().getBasicSamHeader());
		realigned.setReadName(w.list.get(0).getReadHeader());
		realigned.setReadBases(B("AAACCCGGGT"));
		realigned.setCigarString("1S5M4S");
		realigned.setBaseQualities(new byte[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9} );
		realigned.setReadNegativeStrandFlag(true);
		fw.writeRealignment(realigned, 1);
		assertEquals(3, w.list.size());
		assertEquals(3, w.list.size());
		assertEquals("ACCC", w.list.get(1).getReadString());
		assertEquals(0, BreakpointFastqEncoding.getEncodedBreakendOffset(w.list.get(1).getReadHeader()));
		assertEquals("T", w.list.get(2).getReadString());
		assertEquals(9, BreakpointFastqEncoding.getEncodedBreakendOffset(w.list.get(2).getReadHeader()));
		fw.close();
	}
}
