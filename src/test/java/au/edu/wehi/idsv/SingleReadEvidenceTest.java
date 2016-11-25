package au.edu.wehi.idsv;

import static org.junit.Assert.*;

import java.util.List;

import org.junit.Test;

import au.edu.wehi.idsv.picard.InMemoryReferenceSequenceFile;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;

public class SingleReadEvidenceTest extends TestHelper {
	@Test
	public void should_not_create_indel_for_XNX_placeholder() {
		for (SAMRecord r : new SAMRecord[] {
				Read(0, 1, "1X10S"),
				Read(0, 1, "2X10S"),
				Read(0, 1, "1X1N1X10S"),
				Read(0, 1, "1X100N1X10S"),
		}) {
			List<SingleReadEvidence> list = SingleReadEvidence.createEvidence(SES(), r);
			assertEquals(1, list.size());
			assertTrue(list.get(0) instanceof SoftClipEvidence);
		}
	}
	@Test
	public void should_not_calculate_homology_if_untemplated_sequence_exists() {
		SAMRecord realign = Read(0, 40, "10S10M");
		realign.setAttribute("SA", "polyA,10,+,9M1S10S,0,0");
		SplitReadEvidence e = SplitReadEvidence.create(SES(), realign).get(0);
		assertEquals("", e.getHomologySequence());
		assertEquals(0, e.getHomologyAnchoredBaseCount());
	}
	@Test
	public void should_not_calculate_homology_if_inexact() {
		SAMRecord realign = Read(0, 1, "1X10S");
		realign.setAttribute("SA", "polyA,10,+,10M,0,0");
		SplitReadEvidence e = SplitReadEvidence.create(SES(), realign).get(0);
		assertFalse(e.isBreakendExact());
		assertEquals("", e.getHomologySequence());
		assertEquals(0, e.getHomologyAnchoredBaseCount());
	}
	public RealignedBreakpoint test_seq(String originalBreakpointSequence, String cigar, BreakendDirection direction, boolean alignNegativeStrand, String expectedUntemplatedSequence) {
		SAMRecord r = Read(0, 1, cigar);
		r.setReadBases(alignNegativeStrand ? B(SequenceUtil.reverseComplement(originalBreakpointSequence)): B(originalBreakpointSequence));
		r.setReadNegativeStrandFlag(alignNegativeStrand);
		RealignedBreakpoint rbp = RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(), new BreakendSummary(0, direction, 1), B("N"), r);
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
		SAMRecord realign = Read(1, 40, "5S2M");
		SplitReadEvidence e = (SplitReadEvidence)incorporateRealignment(SES(), realign, withSequence("NTACG", Read(1, 15, "5M")));
		assertEquals("TACG", e.getHomologySequence());
		assertEquals(4, e.getHomologyAnchoredBaseCount());
		assertEquals(15, e.getBreakendSummary().start);
		assertEquals(19, e.getBreakendSummary().end);
		assertEquals(36, e.getBreakendSummary().start2);
		assertEquals(40, e.getBreakendSummary().end2);
	}
	@Test
	public void microhomology_match_ff() {
		// 0        1         2         3         4         5         6         7         
		// 1234567890123456789012345678901234567890123456789012345678901234567890123456789
		// CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACC
		//               >>>>>|||        >>>>>          
		//               GCGGGTTGGC      GCCAACCCGC
		SAMRecord realign = Read(2, 31, "5M");
		realign.setReadNegativeStrandFlag(true);
		realign.setReadBases(B("GCCAA"));
		SplitReadEvidence e = (SplitReadEvidence)incorporateRealignment(SES(), realign, withSequence("GCGGG", Read(1, 15, "5M")));
		assertEquals("TTG", e.getHomologySequence());
		assertEquals(0, e.getHomologyAnchoredBaseCount());
		assertEquals(19, e.getBreakpointSummary().start);
		assertEquals(22, e.getBreakpointSummary().end);
		assertEquals(32, e.getBreakpointSummary().start2);
		assertEquals(35, e.getBreakpointSummary().end2);
	}
	@Test
	public void microhomology_match_bb() {
		// 0        1         2         3         4         5         6         7         
		// 1234567890123456789012345678901234567890123456789012345678901234567890123456789
		// CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACC
		//       <<<<<                       <<<<<
		//  TCAGTTCGCA                  AGCGAACTGA
		//       ^                          ^ 1bp micro
		SAMRecord realign = Read(2, 35, "5M");
		realign.setReadBases(B("ACTGA"));
		realign.setReadNegativeStrandFlag(true);
		RealignedBreakpoint rbp = RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(), new BreakendSummary(2, BWD, 7), "TCGCA", realign);
		assertEquals("T", rbp.getHomologySequence());
		assertEquals(1, rbp.getHomologyAnchoredBaseCount());
		assertEquals(7, rbp.getBreakpointSummary().start);
		assertEquals(8, rbp.getBreakpointSummary().end);
		assertEquals(34, rbp.getBreakpointSummary().start2);
		assertEquals(35, rbp.getBreakpointSummary().end2);
	}
	@Test
	public void microhomology_match_bf() {
		// 0        1         2         3         4         5         6         7         
		// 1234567890123456789012345678901234567890123456789012345678901234567890123456789
		// ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA
		//          |>****              <----|                
		SAMRecord realign = Read(1, 10, "2M");
		realign.setReadBases(B("NN"));
		RealignedBreakpoint rbp = RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(), new BreakendSummary(1, BWD, 30), "TACGN", realign);
		assertEquals("TACG", rbp.getHomologySequence());
		assertEquals(4, rbp.getHomologyAnchoredBaseCount());
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
		SoftClipEvidence sce = SoftClipEvidence.create(SES(), FWD, r);
		RealignedBreakpoint rbp = RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(), sce.getBreakendSummary(), r.getReadBases(), realign);
		// breakpoint could be anywhere in the poly A microhomology
		assertEquals(rbp.getHomologySequence().length(), rbp.getBreakpointSummary().end - rbp.getBreakpointSummary().start);
		assertEquals(rbp.getHomologySequence().length(), rbp.getBreakpointSummary().end2 - rbp.getBreakpointSummary().start2);
	}
	@Test
	public void microhomology_bases_should_ignore_case() {
		RealignedBreakpoint rbp = RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(), new BreakendSummary(0, FWD, 10), "aaaa", Read(0, 100, "5M"));
		// microhomology in both directions
		assertEquals(9, rbp.getHomologySequence().length());
	}
	@Test
	public void microhomology_bases_should_match_on_appropriate_strand() {
		SAMRecord realign = Read(0, 100, "5M");
		realign.setReadBases(B("NNNN"));
		realign.setReadNegativeStrandFlag(true);
		RealignedBreakpoint rbp = RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(), new BreakendSummary(0, FWD, 10), "TTTT", realign);
		assertEquals("TTTT", rbp.getHomologySequence());
		realign.setReadNegativeStrandFlag(false);
		rbp = RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(), new BreakendSummary(0, FWD, 10), "AAAA", realign);
		assertEquals("AAAA", rbp.getHomologySequence());
	}
	@Test
	public void should_calculate_imprecise_breakpoint() {
		assertEquals(new BreakpointSummary(0, FWD, 150, 100, 200, 1, BWD, 250, 200, 300),
				RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(),
						new BreakendSummary(0, FWD, 150, 100, 200), "",
						R(1, 300, "5M", "GTNCA", false)).getBreakpointSummary());
		
		assertEquals(new BreakpointSummary(0, FWD, 150, 100, 200, 1, FWD, 354, 304, 404),
				RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(),
						new BreakendSummary(0, FWD, 150, 100, 200), "",
						R(1, 300, "5M", "GTNCA", true)).getBreakpointSummary());
		
		assertEquals(new BreakpointSummary(0, BWD, 150, 100, 200, 1, FWD, 354, 304, 404),
				RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(),
						new BreakendSummary(0, BWD, 150, 100, 200), "",
						R(1, 300, "5M", "GTNCA", false)).getBreakpointSummary());
		
		assertEquals(new BreakpointSummary(0, BWD, 150, 100, 200, 1, BWD, 250, 200, 300),
				RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(),
						new BreakendSummary(0, BWD, 150, 100, 200), "",
						R(1, 300, "5M", "GTNCA", true)).getBreakpointSummary());
	}
	@Test
	public void should_not_calculate_microhomology_for_imprecise_breakpoint() {
		assertEquals(0, RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(), new BreakendSummary(0, FWD, 150, 100, 200), "", R(0, 100, "5M", "AAAAA", true)).getHomologySequence().length());
	}
	public static SAMRecord R(final int referenceIndex, final int alignmentStart, final String cigar, final String bases, final boolean negativeStrand) {
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
		assertEquals("GT", RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(), new BreakendSummary(0, FWD, 150, 100, 200), "", R(1, 100, "2S3M", "GTNCA", false)).getInsertedSequence()); 
		assertEquals("GT", RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(), new BreakendSummary(0, FWD, 150, 100, 200), "", R(1, 100, "3M2S", SequenceUtil.reverseComplement("GTNCA"), true)).getInsertedSequence());
		assertEquals("CA", RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(), new BreakendSummary(0, BWD, 150, 100, 200), "", R(1, 100, "3M2S", "GTNCA", false)).getInsertedSequence()); 
		assertEquals("CA", RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(), new BreakendSummary(0, BWD, 150, 100, 200), "", R(1, 100, "2S3M", SequenceUtil.reverseComplement("GTNCA"), true)).getInsertedSequence());
	}
	@Test(expected=IllegalArgumentException.class)
	public void breakpoint_interval_anchor_sequence_should_be_sane() {
		// Can't have inexact breakpoint with anchored sequence
		// (RealignedBreakpoint does the microhomology calculation)
		RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(), new BreakendSummary(0, FWD, 150, 100, 200), "T", Read(0, 40, "2M"));
	}
	@Test
	public void should_consider_microhomology_for_both_breakends() {
		// 12345678901234567890
		// ACGTACGTACGTACGTACGT
		// NCG    TAN
		//    ^^--
		assertEquals(new BreakpointSummary(1, FWD, 3, 1, 5, 1, BWD, 8, 6, 10), RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(), new BreakendSummary(1, FWD, 3, 3, 3), "NCG", R(1, 8, "3M", "TAN", false)).getBreakpointSummary());
	}
	@Test
	public void should_calc_microhomology() {
		InMemoryReferenceSequenceFile ref = new InMemoryReferenceSequenceFile(new String[] { "0", "1" }, new byte[][] {
				B("AAAAAAAA"), B("AAAAAAAA")});
		//       >>>>            <<<<
		RealignedBreakpoint bp = RealignedBreakpoint.create(ref, ref.getSequenceDictionary(), new BreakendSummary(0, FWD, 4), "AAAA", R(1, 5, "4M", "AAAA", false));
		assertEquals("AAAAAAAA", bp.getHomologySequence());
		assertEquals(4, bp.getHomologyAnchoredBaseCount());
		assertEquals(new BreakpointSummary(0, FWD, 4, 1, 8, 1, BWD, 4, 1, 8), bp.getBreakpointSummary());
		
		ref = new InMemoryReferenceSequenceFile(new String[] { "0", "1" }, new byte[][] {
				B("AAAAAAAA"), B("AAAAAAAA")});
		//       >>>            <<<<<
		bp = RealignedBreakpoint.create(ref, ref.getSequenceDictionary(), new BreakendSummary(0, FWD, 3), "AAA", R(1, 4, "5M", "AAAAA", false));
		assertEquals("AAAAAAAA", bp.getHomologySequence());
		assertEquals(3, bp.getHomologyAnchoredBaseCount());
		assertEquals(new BreakpointSummary(0, FWD, 3, 1, 8, 1, BWD, 3, 1, 8), bp.getBreakpointSummary());
		
		ref = new InMemoryReferenceSequenceFile(new String[] { "0", "1" }, new byte[][] {
				B("ATTAGGCG"), B("CGCCTAAT")});
		//       >>>>        >>>>
		bp = RealignedBreakpoint.create(ref, ref.getSequenceDictionary(), new BreakendSummary(0, FWD, 4), "ATTA", R(1, 1, "4M", "CGCC", true));
		assertEquals(new BreakpointSummary(0, FWD, 4, 1, 8, 1, FWD, 4, 1, 8), bp.getBreakpointSummary());
		assertEquals("ATTAGGCG", bp.getHomologySequence());
		assertEquals(4, bp.getHomologyAnchoredBaseCount());
		
		ref = new InMemoryReferenceSequenceFile(new String[] { "0", "1" }, new byte[][] {
				B("NNNNAAAA"), B("NNAAAA")});
		//       >>>>          <<<<
		bp = RealignedBreakpoint.create(ref, ref.getSequenceDictionary(), new BreakendSummary(0, FWD, 4), "NNNN", R(1, 3, "4M", "AAAA", false));
		assertEquals("AAAA", bp.getHomologySequence());
		assertEquals(0, bp.getHomologyAnchoredBaseCount());
		
		ref = new InMemoryReferenceSequenceFile(new String[] { "0", "1" }, new byte[][] {
				B("TAAANNNN"), B("NNTTTAN")});
		//           <<<<      <<<<    
		bp = RealignedBreakpoint.create(ref, ref.getSequenceDictionary(), new BreakendSummary(0, BWD, 5), "NNNN", R(1, 3, "4M", "TTTA", true));
		assertEquals("TAAA", bp.getHomologySequence());
		assertEquals(0, bp.getHomologyAnchoredBaseCount());
		
		ref = new InMemoryReferenceSequenceFile(new String[] { "0", "1" }, new byte[][] {
				B("TAAANNNN"), B("NNTAAAN")});
		//           <<<<      >>>>    
		bp = RealignedBreakpoint.create(ref, ref.getSequenceDictionary(), new BreakendSummary(0, BWD, 5), "NNNN", R(1, 3, "4M", "TAAA", false));
		assertEquals("TAAA", bp.getHomologySequence());
		assertEquals(0, bp.getHomologyAnchoredBaseCount());
	}
	@Test
	public void should_call_deletion_microhomology() {
		//          1         2         3         4
		// 1234567890123456789012345678901234567890123456789
		// CATTAATCGCAATAAAATGTTCAAAACGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAA
		// ***************         ***************
		String contig = "CATTAATCGCAAGAAAAGGTTGAAAACGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAA";
		InMemoryReferenceSequenceFile ref = new InMemoryReferenceSequenceFile(new String[] { "Contig" }, new byte[][] { B(contig) });
		RealignedBreakpoint bp = RealignedBreakpoint.create(ref, ref.getSequenceDictionary(), new BreakendSummary(0, FWD, 15), "CATTAATCGCAATAA", R(0, 25, "25M", "AACGACGCCAAGTCA", false));
		assertEquals("AAAA", bp.getHomologySequence());
	}
	@Test
	public void inserted_sequence_should_be_local_breakend_strand() {
		String contig = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
		InMemoryReferenceSequenceFile ref = new InMemoryReferenceSequenceFile(new String[] { "Contig" }, new byte[][] { B(contig) });
		RealignedBreakpoint rbp = RealignedBreakpoint.create(ref, ref.getSequenceDictionary(), new BreakendSummary(0, FWD, 15), "N", onNegative(R(0, 1, "5M5S", "AAAAGGTTA", false))[0]);
		assertEquals("TAACC", rbp.getInsertedSequence());
	}
}
