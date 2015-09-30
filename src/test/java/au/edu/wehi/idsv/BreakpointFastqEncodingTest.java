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
		public int getLocalMaxBaseQual() {
			return 4;
		}
		@Override
		public int getLocalTotalBaseQual() {
			return 5;
		}
		@Override
		public boolean isBreakendExact() {
			return false;
		}
		@Override
		public float getBreakendQual() {
			return 6;
		}
		@Override
		public byte[] getAnchorSequence() {
			// TODO Auto-generated method stub
			return B("A");
		}
		@Override
		public byte[] getAnchorQuality() {
			// TODO Auto-generated method stub
			return B("A");
		}
	}
	@Test
	public void should_round_trip_ID() {
		FastqRecord fq = BreakpointFastqEncoding.getRealignmentFastq(new StubDirectedBreakend());
		assertEquals("EvidenceID", BreakpointFastqEncoding.getEncodedID(fq.getReadHeader()));
	}
	@Test
	public void should_round_trip_StartPosition() {
		FastqRecord fq = BreakpointFastqEncoding.getRealignmentFastq(new StubDirectedBreakend());
		assertEquals(123, BreakpointFastqEncoding.getEncodedStartPosition(fq.getReadHeader()));
	}
	@Test
	public void should_round_trip_ReferenceIndex() {
		FastqRecord fq = BreakpointFastqEncoding.getRealignmentFastq(new StubDirectedBreakend());
		assertEquals(1, BreakpointFastqEncoding.getEncodedReferenceIndex(fq.getReadHeader()));
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
	@Test
	public void soft_clip_should_be_sorted_by_sam_coordinate_sort_order() {
		FastqRecord fq = BreakpointFastqEncoding.getRealignmentFastq(SCE(FWD, Read(1, 2, "5M5S")));
		assertEquals(1, BreakpointFastqEncoding.getEncodedReferenceIndex(fq.getReadHeader()));
		assertEquals(2, BreakpointFastqEncoding.getEncodedStartPosition(fq.getReadHeader()));
	}
	@Test
	public void sam_assembly_should_be_sorted_by_sam_coordinate_sort_order() {
		SAMRecordAssemblyEvidence ass =  AssemblyFactory.createAnchoredBreakend(getContext(), AES(), BWD, null, 1, 4, 3, B("AAAAAA"), B("AAAAAA"));
		FastqRecord fq = BreakpointFastqEncoding.getRealignmentFastq(ass);
		assertEquals(1, BreakpointFastqEncoding.getEncodedReferenceIndex(fq.getReadHeader()));
		assertEquals(ass.getSAMRecord().getAlignmentStart(), BreakpointFastqEncoding.getEncodedStartPosition(fq.getReadHeader()));
	}
}
