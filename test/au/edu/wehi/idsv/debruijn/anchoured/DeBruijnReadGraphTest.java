package au.edu.wehi.idsv.debruijn.anchoured;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import htsjdk.samtools.SAMRecord;

import org.junit.Test;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.SoftClipEvidence;
import au.edu.wehi.idsv.debruijn.anchored.DeBruijnReadGraph;
import au.edu.wehi.idsv.sam.AnomolousReadAssembly;
import au.edu.wehi.idsv.TestHelper;

public class DeBruijnReadGraphTest extends TestHelper {
	private static SAMRecord R(String read) {
		return R(null, read, null, false, true);
	}
	private static SAMRecord R(String readName, String read) {
		return R(readName, read, null, false, true);
	}
	private static SAMRecord R(String readName, String read, byte[] qual, boolean mappedNegativeStrand, boolean mateNegativeStrand) {
		SAMRecord record = new SAMRecord(getHeader());
		if (qual == null) {
			qual = new byte[read.length()];
			for (int i = 0; i < qual.length; i++) qual[i] = 1;
		}
		record.setReadBases(B(read));
		record.setBaseQualities(qual);
		record.setReadPairedFlag(true);
		record.setReadNegativeStrandFlag(mappedNegativeStrand);
		record.setReadUnmappedFlag(false);
		record.setMateUnmappedFlag(false);
		record.setMateNegativeStrandFlag(mateNegativeStrand);
		if (readName == null) {
			readName = String.format("%s-%s-%s%s", read, qual, mappedNegativeStrand, mateNegativeStrand);
		}
		record.setReadName(readName);
		record.setReferenceIndex(0);
		record.setAlignmentStart(1);
		record.setCigarString(String.format("%dM1S", read.length() - 1));
		return record;
	}
	private SAMRecord inferLocal(DeBruijnReadGraph ass, SAMRecord remote) {
		SAMRecord local = Read(0, 1, "1M");
		local.setReadName(remote.getReadName());
		local.setReadPairedFlag(true);
		local.setReadNegativeStrandFlag(ass.getDirection() == BWD);
		return local;
	}
	private void addRead(DeBruijnReadGraph ass, SAMRecord r, boolean sc) {
		if (sc) {
			ass.addEvidence(SCE(ass.getDirection(), r));
		} else {
			ass.addEvidence(NRRP(inferLocal(ass, r), r));
		}
	}
	private void removeRead(DeBruijnReadGraph ass, SAMRecord r, boolean sc) {
		if (sc) {
			ass.removeEvidence(SCE(ass.getDirection(), r));
		} else {
			ass.removeEvidence(NRRP(inferLocal(ass, r), r));
		}
	}
	@Test
	public void should_assemble_single_read() {
		DeBruijnReadGraph ass = new DeBruijnReadGraph(4, BreakendDirection.Forward);
		addRead(ass, R("AAAACGTC"), true);
		SAMRecord result = ass.assembleVariant();
		assertEquals("AAAACGTC", S(result.getReadBases()));
	}
	@Test
	public void should_assemble_positive_strand_consensus() {
		DeBruijnReadGraph ass = new DeBruijnReadGraph(4, BreakendDirection.Forward);
		addRead(ass, R(null, "AAAACGTC", null, true, true), true);
		SAMRecord result = ass.assembleVariant();
		assertEquals("AAAACGTC", S(result.getReadBases()));
	}
	@Test
	public void should_assemble_unanchored_reads() {
		DeBruijnReadGraph ass = new DeBruijnReadGraph(3, BreakendDirection.Forward);
		addRead(ass, withSequence("CTAAA", Read(0, 1, "4M1S"))[0], true);
		addRead(ass, R(null, "AAAGT", null, false, true), false);
		SAMRecord result = ass.assembleVariant();
		assertEquals("CTAAAGT", S(result.getReadBases()));
		assertEquals("4M3S", result.getCigarString());
	}
	@Test(expected = RuntimeException.class)  
	public void unanchored_reads_should_require_mapped_mate() {
		DeBruijnReadGraph ass = new DeBruijnReadGraph(3, BreakendDirection.Forward);
		addRead(ass, R("CTAAA"), true);
		SAMRecord unanchored = R(null, "AAAGT", null, false, true);
		unanchored.setMateUnmappedFlag(true);
		addRead(ass, unanchored, false);
		SAMRecord result = ass.assembleVariant();
		assertEquals("CTAAAGT", S(result.getReadBases()));
		assertEquals("5M2S", result.getCigarString());
	}
	@Test
	public void should_assemble_unanchored_reads_in_FR_orientation() {
		assertExpected(BreakendDirection.Forward, "CTAAAGT", "4M3S", "CTAAA", "AAAGT", true, false);
		assertExpected(BreakendDirection.Forward, "CTAAAGT", "4M3S", "CTAAA", "AAAGT", false, true);
		assertExpected(BreakendDirection.Forward, "CTAAAGT", "4M3S", "CTAAA", "ACTTT", true, true);
		assertExpected(BreakendDirection.Forward, "CTAAAGT", "4M3S", "CTAAA", "ACTTT", false, false);
		
		assertExpected(BreakendDirection.Backward, "GTAAACT", "3S4M", "AAACT", "GTAAA", true, false);
		assertExpected(BreakendDirection.Backward, "GTAAACT", "3S4M", "AAACT", "GTAAA", false, true);
		assertExpected(BreakendDirection.Backward, "GTAAACT", "3S4M", "AAACT", "TTTAC", true, true);
		assertExpected(BreakendDirection.Backward, "GTAAACT", "3S4M", "AAACT", "TTTAC", false, false);
	}
	private void assertExpected(BreakendDirection direction, String expectedSeq, String expectedCigar , String anchorSeq, String unanchorSeq, boolean mappedNegativeStrand, boolean mateNegativeStrand) {
		SAMRecord anchor = R(anchorSeq);
		if (direction == BWD) anchor.setCigarString("1S4M");
		
		// Assembly should not depend on whether the read is mapped or not 
		DeBruijnReadGraph ass = new DeBruijnReadGraph(3, direction);
		addRead(ass, anchor, true);
		SAMRecord unanchored = R(null, unanchorSeq, null, mappedNegativeStrand, mateNegativeStrand);
		unanchored.setReadUnmappedFlag(true);
		addRead(ass, unanchored, false);
		SAMRecord result = ass.assembleVariant();
		assertEquals(expectedSeq, S(result.getReadBases()));
		assertEquals(expectedCigar, result.getCigarString());
		
		ass = new DeBruijnReadGraph(3, direction);
		addRead(ass, anchor, true);
		unanchored = R(null, unanchorSeq, null, mappedNegativeStrand, mateNegativeStrand);
		unanchored.setReadUnmappedFlag(false);
		addRead(ass, unanchored, false);
		result = ass.assembleVariant();
		assertEquals(expectedSeq, S(result.getReadBases()));
		assertEquals(expectedCigar, result.getCigarString());
	}
	@Test
	public void should_assemble_soft_clipped_read() {
		DeBruijnReadGraph ass = new DeBruijnReadGraph(4, BreakendDirection.Forward);
		SAMRecord sc = R("AAAACGTC");
		sc.setCigarString("4M4S");
		addRead(ass, sc, true);
		SAMRecord result = ass.assembleVariant();
		assertEquals("AAAACGTC", S(result.getReadBases()));
		assertEquals("4M4S", result.getCigarString());
	}
	@Test
	public void should_greedy_traverse_highest_weight_path() {
		DeBruijnReadGraph ass = new DeBruijnReadGraph(2, BreakendDirection.Forward);
		SAMRecord sc = R(null, "ACGTACTGAG", new byte[] { 1,2,3,4,5,6,7,8,9,10}, false, true);
		sc.setCigarString("4M6S");
		addRead(ass, sc, true);
		SAMRecord result = ass.assembleVariant();
		assertEquals("TACTGAGT", S(result.getReadBases()));
	}
	@Test
	public void should_assemble_backward_breakpoint() {
		DeBruijnReadGraph ass = new DeBruijnReadGraph(4, BreakendDirection.Backward);
		SAMRecord sc = R("AAAACGTC");
		sc.setCigarString("4S4M");
		addRead(ass, sc, true);
		SAMRecord result = ass.assembleVariant();
		assertEquals("AAAACGTC", S(result.getReadBases()));
		assertEquals("4S4M", result.getCigarString());
	}
	@Test
	public void remove_should_exclude_from_graph() {
		DeBruijnReadGraph ass = new DeBruijnReadGraph(4, BreakendDirection.Forward);
		SAMRecord sc = R("AAAACGTC");
		sc.setCigarString("4M4S");
		addRead(ass, sc, true);
		addRead(ass, R("r1", "AAAATTTT"), false);
		addRead(ass, R("r2", "AAAATTTT"), false);
		assertEquals("AAAATTTT", S(ass.assembleVariant().getReadBases()));
		removeRead(ass, R("r1", "AAAATTTT"), false);
		removeRead(ass, R("r2", "AAAATTTT"), false);
		SAMRecord result = ass.assembleVariant();
		assertEquals("AAAACGTC", S(result.getReadBases()));
	}
	@Test
	public void should_use_offset_kmer_if_softclip_longer_than_k() {
		DeBruijnReadGraph ass = new DeBruijnReadGraph(4, BreakendDirection.Forward);
		SAMRecord sc = R("AAAACGTC");
		sc.setCigarString("2M6S");
		addRead(ass, sc, true);
		SAMRecord result = ass.assembleVariant();
		assertEquals("AAAACGTC", S(result.getReadBases()));
		assertEquals("2M6S", result.getCigarString());
	}
	@Test
	public void assembly_base_quality_should_be_sum_of_min_read_qualities() {
		DeBruijnReadGraph ass = new DeBruijnReadGraph(3, BreakendDirection.Forward);
		addRead(ass, R(null, "ACGTA", new byte[] { 1,2,3,4,5 }, false, true), true);
		addRead(ass, R(null, "ACGTA", new byte[] { 3,4,5,6,7 }, false, true), true);
		SAMRecord result = ass.assembleVariant();
		// pad out read qualities
		assertArrayEquals(new byte[] { 4,6,8,8,8 }, result.getBaseQualities());
	}
	@Test
	public void assembly_base_quality_should_pad_to_match_read_length() {
		DeBruijnReadGraph ass = new DeBruijnReadGraph(3, BreakendDirection.Forward);
		addRead(ass, R(null, "ACGTA", new byte[] { 1,2,3,4,5 }, false, true), true);
		addRead(ass, R(null, "ACGTA", new byte[] { 3,4,5,6,7 }, false, true), true);
		SAMRecord result = ass.assembleVariant();
		// pad out read qualities
		assertEquals(result.getReadBases().length, result.getBaseQualities().length);
	}
	@Test
	public void read_count_should_be_number_of_reads__with_at_least_one_kmer_on_path() {
		DeBruijnReadGraph ass = new DeBruijnReadGraph(3, BreakendDirection.Forward);
		addRead(ass, R(null, "ACGTT", new byte[] { 1,2,3,4,5 }, false, true), true);
		addRead(ass, R(null, "ACGTA", new byte[] { 3,4,5,6,7 }, false, true), true);
		addRead(ass, R(null, "AAAAA", new byte[] { 1,1,1,1,1 }, false, true), true);
		addRead(ass, R(null, "ACGAA", new byte[] { 1,1,1,1,1 }, false, true), true);
		AnomolousReadAssembly result = ass.assembleVariant();
		
		assertEquals(3, (int)result.getReadCount());
	}
	// hack so we don't have to redo the test case as the test was created
	// with 0 length soft clips
	public SoftClipEvidence HackSCE(String seq, byte[] qual) {
		SAMRecord read = withQual(qual, withSequence(seq, Read(0, 1, "4M1S")))[0];
		SoftClipEvidence e = new SoftClipEvidence(getContext(), FWD, read);
		read.setCigarString("5M");
		return e;
	}
	@Test
	public void read_base_count_should_be_number_of_read_bases_on_returned_path() {
		DeBruijnReadGraph ass = new DeBruijnReadGraph(3, BreakendDirection.Forward);
		ass.addEvidence(HackSCE("ACGTT", new byte[] { 1,2,3,4,5 })); // 4
		ass.addEvidence(HackSCE("ACGTA", new byte[] { 3,4,5,6,7 })); // 5
		ass.addEvidence(HackSCE("AAAAA", new byte[] { 1,1,1,1,1 })); // 0 since no kmer on path
		ass.addEvidence(HackSCE("ACGAA", new byte[] { 1,1,1,1,1 })); // 3
		addRead(ass, R(null, "TTGTA", new byte[] { 1,1,1,1,1 }, false, true), false); // 3
		addRead(ass, R(null, "GTACG", new byte[] { 1,1,1,1,1 }, false, true), false); // 4+3 since kmers are disconnected 
		AnomolousReadAssembly result = ass.assembleVariant();
		
		assertEquals("test assumes this contruction - did I get the test case wrong?", "TACGTA", S(result.getReadBases()));
		assertEquals(4 + 5 + 3 + 3 + 4+3, (int)result.getReadBaseCount());
	}
}