package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;

import org.junit.Test;


public class RealignedBreakpointTest extends TestHelper {
	@Test(expected=IllegalArgumentException.class)
	public void should_throw_if_realigned_unmapped() {
		AssemblyEvidence dba = AE();
		new RealignedBreakpoint(getContext(), dba.getBreakendSummary(), B("N"), Unmapped(2));
	}
	@Test
	public void should_set_realign_evidence() {
		new RealignedBreakpoint(getContext(), new BreakendSummary(0, FWD, 1, 1), B("N"), withMapq(10, Read(0, 1, "5M"))[0]);
	}
	public RealignedBreakpoint test_seq(String originalBreakpointSequence, String cigar, BreakendDirection direction, boolean alignNegativeStrand, String expectedUntemplatedSequence) {
		SAMRecord r = Read(0, 1, cigar);
		r.setReadBases(alignNegativeStrand ? B(SequenceUtil.reverseComplement(originalBreakpointSequence)): B(originalBreakpointSequence));
		r.setReadNegativeStrandFlag(alignNegativeStrand);
		RealignedBreakpoint rbp = new RealignedBreakpoint(getContext(), new BreakendSummary(0, direction, 1, 1), B("N"), r);
		assertEquals(expectedUntemplatedSequence, rbp.getInsertedSequence());
		return rbp;
	}
	@Test
	public void inserted_sequence_should_be_relative_to_original_realigned_sequence() {
		test_seq("ATTCNNAGC", "4S3M3S", FWD, false, "ATTC");
		test_seq("ATTCNNAGC", "4S3M3S", FWD, true, "ATT");
		test_seq("ATTCNNAGC", "4S3M3S", BWD, false, "AGC");
		test_seq("ATTCNNAGC", "4S3M3S", BWD, true, "NAGC");
	}
	@Test
	public void microhomology_match_fb() {
		// 0        1         2         3         4         5         6         7         
		// 1234567890123456789012345678901234567890123456789012345678901234567890123456789
		// ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA
		//              |---->                ****<| 
		SAMRecord realign = Read(1, 40, "2M");
		RealignedBreakpoint rbp = new RealignedBreakpoint(getContext(), new BreakendSummary(1, FWD, 19, 19), "NTACG", realign);
		assertEquals(4, rbp.getMicroHomologyLength());
		assertEquals(15, rbp.getBreakpointSummary().start);
		assertEquals(19, rbp.getBreakpointSummary().end);
		assertEquals(36, rbp.getBreakpointSummary().start2);
		assertEquals(40, rbp.getBreakpointSummary().end2);
	}
	@Test
	public void microhomology_match_ff() {
		// 0        1         2         3         4         5         6         7         
		// 1234567890123456789012345678901234567890123456789012345678901234567890123456789
		// ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA
		//              |---->                    |>**** 
		SAMRecord realign = Read(1, 40, "2M");
		realign.setReadNegativeStrandFlag(true);
		RealignedBreakpoint rbp = new RealignedBreakpoint(getContext(), new BreakendSummary(1, FWD, 19, 19), "NATGC", realign);
		assertEquals(4, rbp.getMicroHomologyLength());
		assertEquals(15, rbp.getBreakpointSummary().start);
		assertEquals(19, rbp.getBreakpointSummary().end);
		assertEquals(41, rbp.getBreakpointSummary().start2);
		assertEquals(45, rbp.getBreakpointSummary().end2);
	}
	@Test
	public void microhomology_match_bb() {
		// 0        1         2         3         4         5         6         7         
		// 1234567890123456789012345678901234567890123456789012345678901234567890123456789
		// ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA
		//      ****<|                  <----|                
		SAMRecord realign = Read(1, 10, "2M");
		realign.setReadNegativeStrandFlag(true);
		RealignedBreakpoint rbp = new RealignedBreakpoint(getContext(), new BreakendSummary(1, BWD, 30, 30), "aTgCN", realign);
		assertEquals(4, rbp.getMicroHomologyLength());
		assertEquals(30, rbp.getBreakpointSummary().start);
		assertEquals(34, rbp.getBreakpointSummary().end);
		assertEquals(6, rbp.getBreakpointSummary().start2);
		assertEquals(10, rbp.getBreakpointSummary().end2);
	}
	@Test
	public void microhomology_match_bf() {
		// 0        1         2         3         4         5         6         7         
		// 1234567890123456789012345678901234567890123456789012345678901234567890123456789
		// ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA
		//          |>****              <----|                
		SAMRecord realign = Read(1, 10, "2M");
		RealignedBreakpoint rbp = new RealignedBreakpoint(getContext(), new BreakendSummary(1, BWD, 30, 30), "TACGN", realign);
		assertEquals(4, rbp.getMicroHomologyLength());
		assertEquals(30, rbp.getBreakpointSummary().start);
		assertEquals(34, rbp.getBreakpointSummary().end);
		assertEquals(11, rbp.getBreakpointSummary().start2);
		assertEquals(15, rbp.getBreakpointSummary().end2);
	}
	@Test
	public void breakpoint_window_size_should_correspond_to_microhomology_length() {
		SAMRecord r = Read(0, 1, "10M5S");
		r.setReadBases(B("TTTTTAAAAAAAAAT"));
		SAMRecord realign = Read(0, 100, "5M");
		realign.setReadBases(B("AAAAT"));
		SoftClipEvidence sce = new SoftClipEvidence(getContext(), SES(), FWD, r);
		RealignedBreakpoint rbp = new RealignedBreakpoint(getContext(), sce.getBreakendSummary(), r.getReadBases(), realign);
		// breakpoint could be anywhere in the poly A microhomology
		assertEquals(rbp.getMicroHomologyLength(), rbp.getBreakpointSummary().end - rbp.getBreakpointSummary().start);
		assertEquals(rbp.getMicroHomologyLength(), rbp.getBreakpointSummary().end2 - rbp.getBreakpointSummary().start2);
	}
	@Test
	public void microhomology_bases_should_ignore_case() {
		RealignedBreakpoint rbp = new RealignedBreakpoint(getContext(), new BreakendSummary(0, FWD, 10, 10), "aaaa", Read(0, 100, "5M"));
		assertEquals(4, rbp.getMicroHomologyLength());
	}
	@Test
	public void microhomology_bases_should_match_on_appropriate_strand() {
		// Turns out strand doesn't really matter because both anchor and realign bases are
		// +ve strand in the SAMRecord, so the complimenting on negative strand has already been done
		SAMRecord realign = Read(0, 100, "5M");
		realign.setReadBases(B("AAAAA"));
		realign.setReadNegativeStrandFlag(true);
		RealignedBreakpoint rbp = new RealignedBreakpoint(getContext(), new BreakendSummary(0, FWD, 10, 10), "AAAA", realign);
		assertEquals(4, rbp.getMicroHomologyLength());
		realign.setReadNegativeStrandFlag(false);
		rbp = new RealignedBreakpoint(getContext(), new BreakendSummary(0, FWD, 10, 10), "AAAA", realign);
		assertEquals(4, rbp.getMicroHomologyLength());
	}
	@Test
	public void should_calculate_imprecise_breakpoint() {
		assertEquals(new BreakpointSummary(0, FWD, 100, 200, 1, BWD, 200, 300),
				new RealignedBreakpoint(getContext(),
						new BreakendSummary(0, FWD, 100, 200), "",
						R(1, 300, "5M", "GTNCA", false)).getBreakpointSummary());
		
		assertEquals(new BreakpointSummary(0, FWD, 100, 200, 1, FWD, 304, 404),
				new RealignedBreakpoint(getContext(),
						new BreakendSummary(0, FWD, 100, 200), "",
						R(1, 300, "5M", "GTNCA", true)).getBreakpointSummary());
		
		assertEquals(new BreakpointSummary(0, BWD, 100, 200, 1, FWD, 304, 404),
				new RealignedBreakpoint(getContext(),
						new BreakendSummary(0, BWD, 100, 200), "",
						R(1, 300, "5M", "GTNCA", false)).getBreakpointSummary());
		
		assertEquals(new BreakpointSummary(0, BWD, 100, 200, 1, BWD, 200, 300),
				new RealignedBreakpoint(getContext(),
						new BreakendSummary(0, BWD, 100, 200), "",
						R(1, 300, "5M", "GTNCA", true)).getBreakpointSummary());
	}
	@Test
	public void should_not_calculate_microhomology_for_imprecise_breakpoint() {
		assertEquals(0, new RealignedBreakpoint(getContext(), new BreakendSummary(0, FWD, 100, 200), "", R(0, 100, "5M", "AAAAA", true)).getMicroHomologyLength());
	}
	private SAMRecord R(final int referenceIndex, final int alignmentStart, final String cigar, final String bases, final boolean negativeStrand) {
		return new SAMRecord(getContext().getBasicSamHeader()) {{
			setReferenceIndex(referenceIndex);
			setAlignmentStart(alignmentStart);
			setReadBases(B(bases));
			setCigarString(cigar);
			setReadNegativeStrandFlag(negativeStrand);
		}};
	}
	@Test
	public void should_include_untemplated_sequence_for_imprecise_breakpoint() {
		assertEquals("GT", new RealignedBreakpoint(getContext(), new BreakendSummary(0, FWD, 100, 200), "", R(1, 100, "3M2S", "GTNCA", false)).getInsertedSequence()); 
		assertEquals("TG", new RealignedBreakpoint(getContext(), new BreakendSummary(0, FWD, 100, 200), "", R(1, 100, "2S3M", "GTNCA", true)).getInsertedSequence());
		assertEquals("CA", new RealignedBreakpoint(getContext(), new BreakendSummary(0, BWD, 100, 200), "", R(1, 100, "2S3M", "GTNCA", false)).getInsertedSequence()); 
		assertEquals("AC", new RealignedBreakpoint(getContext(), new BreakendSummary(0, BWD, 100, 200), "", R(1, 100, "3M2S", "GTNCA", true)).getInsertedSequence());
	}
	@Test(expected=IllegalArgumentException.class)
	public void breakpoint_interval_anchor_sequence_should_be_sane() {
		// Can't have inexact breakpoint with anchored sequence
		// (RealignedBreakpoint does the microhomology calculation)
		new RealignedBreakpoint(getContext(), new BreakendSummary(0, FWD, 100, 200), "T", Read(0, 40, "2M"));
	}
}
