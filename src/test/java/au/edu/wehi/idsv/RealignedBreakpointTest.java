package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.SequenceUtil;

import java.io.IOException;
import java.util.Arrays;

import org.junit.Test;


public class RealignedBreakpointTest extends TestHelper {
	@Test(expected=IllegalArgumentException.class)
	public void should_throw_if_realigned_unmapped() {
		AssemblyEvidence dba = AE();
		RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(), dba.getBreakendSummary(), B("N"), Unmapped(2));
	}
	@Test
	public void should_set_realign_evidence() {
		RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(), new BreakendSummary(0, FWD, 1, 1), B("N"), withMapq(10, Read(0, 1, "5M"))[0]);
	}
	public RealignedBreakpoint test_seq(String originalBreakpointSequence, String cigar, BreakendDirection direction, boolean alignNegativeStrand, String expectedUntemplatedSequence) {
		SAMRecord r = Read(0, 1, cigar);
		r.setReadBases(alignNegativeStrand ? B(SequenceUtil.reverseComplement(originalBreakpointSequence)): B(originalBreakpointSequence));
		r.setReadNegativeStrandFlag(alignNegativeStrand);
		RealignedBreakpoint rbp = RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(), new BreakendSummary(0, direction, 1, 1), B("N"), r);
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
		RealignedBreakpoint rbp = RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(), new BreakendSummary(1, FWD, 19, 19), "NTACG", realign);
		assertEquals("TACG", rbp.getHomologySequence());
		assertEquals(4, rbp.getHomologyBaseCountIncludedLocally());
		assertEquals(15, rbp.getBreakpointSummary().start);
		assertEquals(19, rbp.getBreakpointSummary().end);
		assertEquals(36, rbp.getBreakpointSummary().start2);
		assertEquals(40, rbp.getBreakpointSummary().end2);
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
		RealignedBreakpoint rbp = RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(), new BreakendSummary(2, FWD, 19, 19), "GCGGG", realign);
		assertEquals("TTG", rbp.getHomologySequence());
		assertEquals(0, rbp.getHomologyBaseCountIncludedLocally());
		assertEquals(19, rbp.getBreakpointSummary().start);
		assertEquals(22, rbp.getBreakpointSummary().end);
		assertEquals(32, rbp.getBreakpointSummary().start2);
		assertEquals(35, rbp.getBreakpointSummary().end2);
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
		RealignedBreakpoint rbp = RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(), new BreakendSummary(2, BWD, 7, 7), "TCGCA", realign);
		assertEquals("T", rbp.getHomologySequence());
		assertEquals(1, rbp.getHomologyBaseCountIncludedLocally());
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
		RealignedBreakpoint rbp = RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(), new BreakendSummary(1, BWD, 30, 30), "TACGN", realign);
		assertEquals("TACG", rbp.getHomologySequence());
		assertEquals(4, rbp.getHomologyBaseCountIncludedLocally());
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
		SoftClipEvidence sce = new SoftClipEvidence(SES(), FWD, r);
		RealignedBreakpoint rbp = RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(), sce.getBreakendSummary(), r.getReadBases(), realign);
		// breakpoint could be anywhere in the poly A microhomology
		assertEquals(rbp.getHomologySequence().length(), rbp.getBreakpointSummary().end - rbp.getBreakpointSummary().start);
		assertEquals(rbp.getHomologySequence().length(), rbp.getBreakpointSummary().end2 - rbp.getBreakpointSummary().start2);
	}
	@Test
	public void microhomology_bases_should_ignore_case() {
		RealignedBreakpoint rbp = RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(), new BreakendSummary(0, FWD, 10, 10), "aaaa", Read(0, 100, "5M"));
		// microhomology in both directions
		assertEquals(9, rbp.getHomologySequence().length());
	}
	@Test
	public void microhomology_bases_should_match_on_appropriate_strand() {
		SAMRecord realign = Read(0, 100, "5M");
		realign.setReadBases(B("NNNN"));
		realign.setReadNegativeStrandFlag(true);
		RealignedBreakpoint rbp = RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(), new BreakendSummary(0, FWD, 10, 10), "TTTT", realign);
		assertEquals("TTTT", rbp.getHomologySequence());
		realign.setReadNegativeStrandFlag(false);
		rbp = RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(), new BreakendSummary(0, FWD, 10, 10), "AAAA", realign);
		assertEquals("AAAA", rbp.getHomologySequence());
	}
	@Test
	public void should_calculate_imprecise_breakpoint() {
		assertEquals(new BreakpointSummary(0, FWD, 100, 200, 1, BWD, 200, 300),
				RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(),
						new BreakendSummary(0, FWD, 100, 200), "",
						R(1, 300, "5M", "GTNCA", false)).getBreakpointSummary());
		
		assertEquals(new BreakpointSummary(0, FWD, 100, 200, 1, FWD, 304, 404),
				RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(),
						new BreakendSummary(0, FWD, 100, 200), "",
						R(1, 300, "5M", "GTNCA", true)).getBreakpointSummary());
		
		assertEquals(new BreakpointSummary(0, BWD, 100, 200, 1, FWD, 304, 404),
				RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(),
						new BreakendSummary(0, BWD, 100, 200), "",
						R(1, 300, "5M", "GTNCA", false)).getBreakpointSummary());
		
		assertEquals(new BreakpointSummary(0, BWD, 100, 200, 1, BWD, 200, 300),
				RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(),
						new BreakendSummary(0, BWD, 100, 200), "",
						R(1, 300, "5M", "GTNCA", true)).getBreakpointSummary());
	}
	@Test
	public void should_not_calculate_microhomology_for_imprecise_breakpoint() {
		assertEquals(0, RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(), new BreakendSummary(0, FWD, 100, 200), "", R(0, 100, "5M", "AAAAA", true)).getHomologySequence().length());
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
		assertEquals("GT", RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(), new BreakendSummary(0, FWD, 100, 200), "", R(1, 100, "2S3M", "GTNCA", false)).getInsertedSequence()); 
		assertEquals("GT", RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(), new BreakendSummary(0, FWD, 100, 200), "", R(1, 100, "3M2S", SequenceUtil.reverseComplement("GTNCA"), true)).getInsertedSequence());
		assertEquals("CA", RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(), new BreakendSummary(0, BWD, 100, 200), "", R(1, 100, "3M2S", "GTNCA", false)).getInsertedSequence()); 
		assertEquals("CA", RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(), new BreakendSummary(0, BWD, 100, 200), "", R(1, 100, "2S3M", SequenceUtil.reverseComplement("GTNCA"), true)).getInsertedSequence());
	}
	@Test(expected=IllegalArgumentException.class)
	public void breakpoint_interval_anchor_sequence_should_be_sane() {
		// Can't have inexact breakpoint with anchored sequence
		// (RealignedBreakpoint does the microhomology calculation)
		RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(), new BreakendSummary(0, FWD, 100, 200), "T", Read(0, 40, "2M"));
	}
	@Test
	public void should_consider_microhomology_for_both_breakends() {
		// 12345678901234567890
		// ACGTACGTACGTACGTACGT
		// NCG    TAN
		//    ^^--
		assertEquals(new BreakpointSummary(1, FWD, 1, 5, 1, BWD, 6, 10), RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(), new BreakendSummary(1, FWD, 3, 3), "NCG", R(1, 8, "3M", "TAN", false)).getBreakpointSummary());
	}
	@Test
	public void should_calc_microhomology() {
		MockReference ref = new MockReference(
				"AAAAAAAA", "AAAAAAAA");
		//       >>>>            <<<<
		RealignedBreakpoint bp = RealignedBreakpoint.create(ref, ref.getSequenceDictionary(), new BreakendSummary(0, FWD, 4, 4), "AAAA", R(1, 5, "4M", "AAAA", false));
		assertEquals("AAAAAAAA", bp.getHomologySequence());
		assertEquals(4, bp.getHomologyBaseCountIncludedLocally());
		assertEquals(new BreakpointSummary(0, FWD, 1, 8, 1, BWD, 1, 8), bp.getBreakpointSummary());
		
		ref = new MockReference(
				"AAAAAAAA", "AAAAAAAA");
		//       >>>            <<<<<
		bp = RealignedBreakpoint.create(ref, ref.getSequenceDictionary(), new BreakendSummary(0, FWD, 3, 3), "AAA", R(1, 4, "5M", "AAAAA", false));
		assertEquals("AAAAAAAA", bp.getHomologySequence());
		assertEquals(3, bp.getHomologyBaseCountIncludedLocally());
		assertEquals(new BreakpointSummary(0, FWD, 1, 8, 1, BWD, 1, 8), bp.getBreakpointSummary());
		
		ref = new MockReference(
				"ATTAGGCG", "CGCCTAAT");
		//       >>>>        >>>>
		bp = RealignedBreakpoint.create(ref, ref.getSequenceDictionary(), new BreakendSummary(0, FWD, 4, 4), "ATTA", R(1, 1, "4M", "CGCC", true));
		assertEquals(new BreakpointSummary(0, FWD, 1, 8, 1, FWD, 1, 8), bp.getBreakpointSummary());
		assertEquals("ATTAGGCG", bp.getHomologySequence());
		assertEquals(4, bp.getHomologyBaseCountIncludedLocally());
		
		ref = new MockReference(
				"NNNNAAAA", "NNAAAA");
		//       >>>>          <<<<
		bp = RealignedBreakpoint.create(ref, ref.getSequenceDictionary(), new BreakendSummary(0, FWD, 4, 4), "NNNN", R(1, 3, "4M", "AAAA", false));
		assertEquals("AAAA", bp.getHomologySequence());
		assertEquals(0, bp.getHomologyBaseCountIncludedLocally());
		
		ref = new MockReference(
				"TAAANNNN", "NNTTTAN");
		//           <<<<      <<<<    
		bp = RealignedBreakpoint.create(ref, ref.getSequenceDictionary(), new BreakendSummary(0, BWD, 5, 5), "NNNN", R(1, 3, "4M", "TTTA", true));
		assertEquals("TAAA", bp.getHomologySequence());
		assertEquals(0, bp.getHomologyBaseCountIncludedLocally());
		
		ref = new MockReference(
				"TAAANNNN", "NNTAAAN");
		//           <<<<      >>>>    
		bp = RealignedBreakpoint.create(ref, ref.getSequenceDictionary(), new BreakendSummary(0, BWD, 5, 5), "NNNN", R(1, 3, "4M", "TAAA", false));
		assertEquals("TAAA", bp.getHomologySequence());
		assertEquals(0, bp.getHomologyBaseCountIncludedLocally());
	}
	public class MockReference implements ReferenceSequenceFile {
		private String[] seq;
		private SAMSequenceDictionary dict;
		public MockReference(String... seq) {
			this.seq = seq;
			this.dict = new SAMSequenceDictionary();
			int i = 0;
			for (String s : seq) {
				dict.addSequence(new SAMSequenceRecord(Integer.toString(i++), s.length()));
			}
		}
		@Override
		public SAMSequenceDictionary getSequenceDictionary() {
			return dict;
		}
		@Override
		public ReferenceSequence nextSequence() {
			// TODO Auto-generated method stub
			return null;
		}

		@Override
		public void reset() {
		}

		@Override
		public boolean isIndexed() {
			return true;
		}

		@Override
		public ReferenceSequence getSequence(String contig) {
			return getSubsequenceAt(contig, 1, getSequenceDictionary().getSequence(contig).getSequenceLength());
		}

		@Override
		public ReferenceSequence getSubsequenceAt(String contig, long start, long stop) {
			SAMSequenceRecord s = getSequenceDictionary().getSequence(contig);
			return new ReferenceSequence(s.getSequenceName(), s.getSequenceIndex(), Arrays.copyOfRange(B(seq[s.getSequenceIndex()]), (int)start - 1, (int)stop));
		}

		@Override
		public void close() throws IOException {
			// TODO Auto-generated method stub
			
		}
	}
}
