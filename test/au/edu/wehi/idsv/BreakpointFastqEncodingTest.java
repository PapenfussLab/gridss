package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import htsjdk.samtools.fastq.FastqRecord;

import org.junit.Test;

public class BreakpointFastqEncodingTest extends TestHelper {
	public class StubDirectedBreakend implements DirectedEvidence {
		@Override
		public BreakendSummary getBreakendSummary() {
			return new BreakendSummary(1, FWD, 123, 456);
		}
		@Override
		public String getEvidenceID() {
			return "EvidenceID";
		}
		@Override
		public byte[] getBreakendSequence() {
			return B("GTAC");
		}
		@Override
		public byte[] getBreakendQuality() {
			return new byte[] { 1, 2, 3, 4 };
		}
		@Override
		public EvidenceSource getEvidenceSource() {
			return null;
		}
		@Override
		public int getLocalMapq() {
			return 1;
		}
		@Override
		public int getLocalBaseLength() {
			return 2;
		}
		@Override
		public int getLocalBaseCount() {
			return 3;
		}
		@Override
		public int getLocalMaxBaseQual() {
			return 4;
		}
		@Override
		public int getLocalTotalBaseQual() {
			return 5;
		}
	}
	@Test
	public void should_round_trip_ID() {
		FastqRecord fq = BreakpointFastqEncoding.getRealignmentFastq(new StubDirectedBreakend());
		assertEquals("EvidenceID", BreakpointFastqEncoding.getID(fq.getReadHeader()));
	}
	@Test
	public void should_round_trip_StartPosition() {
		FastqRecord fq = BreakpointFastqEncoding.getRealignmentFastq(new StubDirectedBreakend());
		assertEquals(123, BreakpointFastqEncoding.getStartPosition(fq.getReadHeader()));
	}
	@Test
	public void should_round_trip_ReferenceIndex() {
		FastqRecord fq = BreakpointFastqEncoding.getRealignmentFastq(new StubDirectedBreakend());
		assertEquals(1, BreakpointFastqEncoding.getReferenceIndex(fq.getReadHeader()));
	}
	@Test
	public void should_use_BreakpointSequence() {
		FastqRecord fq = BreakpointFastqEncoding.getRealignmentFastq(new StubDirectedBreakend());
		assertEquals("GTAC", fq.getReadString());
	}
	@Test
	public void should_use_BreakpointQuality() {
		FastqRecord fq = BreakpointFastqEncoding.getRealignmentFastq(new StubDirectedBreakend());
		assertEquals(S(new byte[] { 34, 35, 36, 37}), fq.getBaseQualityString());
	}
}
