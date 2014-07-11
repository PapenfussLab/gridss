package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import htsjdk.samtools.fastq.FastqRecord;

import org.junit.Test;

public class BreakpointFastqEncodingTest extends TestHelper {
	public class StubDirectedBreakpoint implements DirectedBreakpoint {
		@Override
		public BreakendSummary getBreakendSummary() {
			return new BreakendSummary(1, BreakendDirection.Forward, 123, 456, null);
		}
		@Override
		public String getEvidenceID() {
			return "EvidenceID";
		}
		@Override
		public byte[] getBreakpointSequence() {
			return B("GTAC");
		}
		@Override
		public byte[] getBreakpointQuality() {
			return new byte[] { 1, 2, 3, 4 };
		}
	}
	@Test
	public void should_round_trip_ID() {
		FastqRecord fq = BreakpointFastqEncoding.getRealignmentFastq(new StubDirectedBreakpoint());
		assertEquals("EvidenceID", BreakpointFastqEncoding.getID(fq.getReadHeader()));
	}
	@Test
	public void should_round_trip_StartPosition() {
		FastqRecord fq = BreakpointFastqEncoding.getRealignmentFastq(new StubDirectedBreakpoint());
		assertEquals(123, BreakpointFastqEncoding.getStartPosition(fq.getReadHeader()));
	}
	@Test
	public void should_round_trip_ReferenceIndex() {
		FastqRecord fq = BreakpointFastqEncoding.getRealignmentFastq(new StubDirectedBreakpoint());
		assertEquals(1, BreakpointFastqEncoding.getReferenceIndex(fq.getReadHeader()));
	}
	@Test
	public void should_use_BreakpointSequence() {
		FastqRecord fq = BreakpointFastqEncoding.getRealignmentFastq(new StubDirectedBreakpoint());
		assertEquals("GTAC", fq.getReadString());
	}
	@Test
	public void should_use_BreakpointQuality() {
		FastqRecord fq = BreakpointFastqEncoding.getRealignmentFastq(new StubDirectedBreakpoint());
		assertEquals(S(new byte[] { 34, 35, 36, 37}), fq.getBaseQualityString());
	}
}
