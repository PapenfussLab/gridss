package au.edu.wehi.idsv;

import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.junit.Ignore;
import org.junit.Test;

import static org.junit.Assert.*;

public class VariantContextDirectedEvidenceTest extends TestHelper {
	@Test
	public void getBreakendSummary_should_handle_f_single_breakend() {
		BreakendSummary loc = new VariantContextDirectedEvidence(getContext(), AES(), minimalVariant().alleles("A", "A.").log10PError(-1).make()).getBreakendSummary();
		assertEquals(0, loc.referenceIndex);
		assertEquals(1, loc.start);
		assertEquals(1, loc.end);
		assertEquals(BreakendDirection.Forward, loc.direction);
	}
	@Test
	public void getBreakendSummary_should_handle_b_single_breakend() {
		BreakendSummary loc = new VariantContextDirectedEvidence(getContext(), AES(), minimalVariant().alleles("A", ".A").log10PError(-1).make()).getBreakendSummary();
		assertEquals(0, loc.referenceIndex);
		assertEquals(1, loc.start);
		assertEquals(1, loc.end);
		assertEquals(BreakendDirection.Backward, loc.direction);
	}
	@Test
	public void getEvidenceID_should_match_vcf_name() {
		VariantContextDirectedEvidence vc = new VariantContextDirectedEvidence(getContext(), AES(), minimalVariant().alleles("A", "AC[polyA:1[").id("test").make());
		assertEquals("test", vc.getEvidenceID());
	}
	@Test()
	public void should_allow_only_single_allele() {
		assertFalse(new VariantContextDirectedEvidence(getContext(), AES(), minimalVariant().alleles("A", "A.", "AA.").id("test").make()).isValid());
	}
	@Test
	public void getBreakpointSequence_should_match_f_single_alt_allele() {
		VariantContextDirectedEvidence vc = new VariantContextDirectedEvidence(getContext(), AES(), minimalVariant().alleles("A", "ACGT.").id("test").make());
		assertEquals("CGT", vc.getBreakpointSequenceString());
		assertArrayEquals(B(vc.getBreakpointSequenceString()), vc.getBreakendSequence());
	}
	@Test
	public void getBreakpointSequence_should_match_b_single_alt_allele() {
		VariantContextDirectedEvidence vc = new VariantContextDirectedEvidence(getContext(), AES(), minimalVariant().alleles("A", ".CGTA").id("test").make());
		assertEquals("CGT", vc.getBreakpointSequenceString());
		assertArrayEquals(B(vc.getBreakpointSequenceString()), vc.getBreakendSequence());
	}
	@Test
	public void getBreakpointSequence_should_return_untemplated_sequence_for_partnered_b_breakend() {
		VariantContextDirectedEvidence vc = new VariantContextDirectedEvidence(getContext(), AES(), minimalVariant().alleles("A", "[polyA[CGTA").id("test").make());
		assertEquals("CGT", vc.getBreakpointSequenceString());
		assertArrayEquals(B(vc.getBreakpointSequenceString()), vc.getBreakendSequence());
	}
	@Test
	public void getBreakpointSequence_should_return_untemplated_sequence_for_partnered_f_breakend() {
		VariantContextDirectedEvidence vc = new VariantContextDirectedEvidence(getContext(), AES(), minimalVariant().alleles("A", "CGTA[polyA[").id("test").make());
		assertEquals("GTA", vc.getBreakpointSequenceString());
		assertArrayEquals(B(vc.getBreakpointSequenceString()), vc.getBreakendSequence());
	}
	@Test
	public void getBreakendSummary_should_handle_partner_contig() {
		VariantContextDirectedEvidence vc = new VariantContextDirectedEvidence(getContext(), AES(), minimalVariant().start(1).stop(2).alleles("AA", "AAAAT[polyACGT:5[").id("test").make());
		BreakpointSummary loc = (BreakpointSummary)vc.getBreakendSummary();
		assertEquals(0, loc.referenceIndex);
		assertEquals(2, loc.start);
		assertEquals(2, loc.end);
		assertEquals(BreakendDirection.Forward, loc.direction);
		assertEquals(1, loc.referenceIndex2);
		assertEquals(5, loc.start2);
		assertEquals(5, loc.end2);
		assertEquals(BreakendDirection.Backward, loc.direction2);
	}
	@Ignore("Named contigs are not yet processed or required by idsv. This test case can be ignored.")
	@Test
	public void should_extend_process_context_sequence_dictionary_when_encountering_named_contig() {
		GenomicProcessingContext pc = getContext();
		int seqCount = pc.getDictionary().size();
		new VariantContextDirectedEvidence(pc, AES(), minimalVariant().chr("<contigNotInReference>").make());
		assertEquals(seqCount + 1, pc.getDictionary().size());
	}
	@Test(expected=IllegalArgumentException.class)
	public void should_not_extend_process_context_sequence_dictionary_when_encountering_named_contig_partner() {
		GenomicProcessingContext pc = getContext();
		int seqCount = pc.getDictionary().size();
		new VariantContextDirectedEvidence(pc, AES(), minimalVariant().start(1).stop(2).alleles("AA", "AAAAT[<contigNotInReference>:1[").id("test").make());
		assertEquals(seqCount + 1, pc.getDictionary().size());
	}
	//@Test
	//public void should_not_extend_process_context_sequence_dictionary_when_encountering_vcf41_placeholder() {
	//	GenomicProcessingContext pc = getContext();
	//	int seqCount = pc.getDictionary().size();
	//	new VariantContextDirectedEvidence(pc, AES(), minimalVariant().start(1).stop(2).alleles("AA", "AAAAT[<UNKNOWN>[").id("test").make());
	//	assertEquals(seqCount, pc.getDictionary().size());
	//}
	@Test
	public void getBreakendSummary_should_handle_ff_partner() {
		VariantContextDirectedEvidence vc = new VariantContextDirectedEvidence(getContext(), AES(), minimalVariant().start(1).stop(2).alleles("AA", "AAAAT]polyACGT:1234]").id("test").make());
		BreakpointSummary loc = (BreakpointSummary)vc.getBreakendSummary();
		assertEquals(0, loc.referenceIndex);
		assertEquals(2, loc.start);
		assertEquals(2, loc.end);
		assertEquals(BreakendDirection.Forward, loc.direction);
		assertEquals(1, loc.referenceIndex2);
		assertEquals(1234, loc.start2);
		assertEquals(1234, loc.end2);
		assertEquals(BreakendDirection.Forward, loc.direction2);
	}
	@Test
	public void getBreakendSummary_should_handle_fb_partner() {
		VariantContextDirectedEvidence vc = new VariantContextDirectedEvidence(getContext(), AES(), minimalVariant().start(1).stop(2).alleles("AA", "AAAAT[polyACGT:1234[").id("test").make());
		BreakpointSummary loc = (BreakpointSummary)vc.getBreakendSummary();
		assertEquals(0, loc.referenceIndex);
		assertEquals(2, loc.start);
		assertEquals(2, loc.end);
		assertEquals(BreakendDirection.Forward, loc.direction);
		assertEquals(1, loc.referenceIndex2);
		assertEquals(1234, loc.start2);
		assertEquals(1234, loc.end2);
		assertEquals(BreakendDirection.Backward, loc.direction2);
	}
	@Test
	public void getBreakendSummary_should_handle_bb_partner() {
		VariantContextDirectedEvidence vc = new VariantContextDirectedEvidence(getContext(), AES(), minimalVariant().start(1).stop(2).alleles("AA", "[polyACGT:1234[AAAAT").id("test").make());
		BreakpointSummary loc = (BreakpointSummary)vc.getBreakendSummary();
		assertEquals(0, loc.referenceIndex);
		assertEquals(1, loc.start);
		assertEquals(1, loc.end);
		assertEquals(BreakendDirection.Backward, loc.direction);
		assertEquals(1, loc.referenceIndex2);
		assertEquals(1234, loc.start2);
		assertEquals(1234, loc.end2);
		assertEquals(BreakendDirection.Backward, loc.direction2);
	}
	@Test
	public void getBreakendSummary_should_handle_bf_partner() {
		VariantContextDirectedEvidence vc = new VariantContextDirectedEvidence(getContext(), AES(), minimalVariant().start(1).stop(2).alleles("AA", "]polyACGT:1234]AAAAT").id("test").make());
		BreakpointSummary loc = (BreakpointSummary)vc.getBreakendSummary();
		assertEquals(0, loc.referenceIndex);
		assertEquals(1, loc.start);
		assertEquals(1, loc.end);
		assertEquals(BreakendDirection.Backward, loc.direction);
		assertEquals(1, loc.referenceIndex2);
		assertEquals(1234, loc.start2);
		assertEquals(1234, loc.end2);
		assertEquals(BreakendDirection.Forward, loc.direction2);
	}
	//@Test
	//public void getBreakendSummary_should_parse_vcf41_compatability_breakpoint_as_breakend() {
	//	// even in 4.2 mode, we should happily parse our backward compatible serialisation as a breakend
	//	VariantContextDirectedEvidence vc = new VariantContextDirectedEvidence(getContext(), AES(), minimalVariant().start(1).stop(1).alleles("A", "A[<UNKNOWN>[").id("test").make());
	//	assertEquals(vc.getBreakendSummary().getClass(), BreakendSummary.class);
	//}
	@Test
	public void matching_breakpoints_should_call_same_evidence() {
		BreakpointSummary s1 = ((VariantContextDirectedBreakpoint)new IdsvVariantContextBuilder(getContext(),
				new VariantContextBuilder()
					.chr("polyA")
					.start(1)
					.stop(1)
					.alleles("A", "A[polyA:10[")
					.make())
			.make()).getBreakendSummary();
		BreakpointSummary s2 = ((VariantContextDirectedBreakpoint)new IdsvVariantContextBuilder(getContext(),
				new VariantContextBuilder()
					.chr("polyA")
					.start(10)
					.stop(10)
					.alleles("A", "]polyA:1]A")
					.make())
			.make()).getBreakendSummary();
		s2 = s2.remoteBreakpoint();
		assertEquals(s1.direction, s2.direction);
		assertEquals(s1.start, s2.start);
		assertEquals(s1.end, s2.end);
		assertEquals(s1.referenceIndex, s2.referenceIndex);
		assertEquals(s1.direction2, s2.direction2);
		assertEquals(s1.start2, s2.start2);
		assertEquals(s1.end2, s2.end2);
		assertEquals(s1.referenceIndex2, s2.referenceIndex2);
	}
	@Test
	public void forward_breakend_should_match_vcf_call() {
		BreakpointSummary s1 = ((VariantContextDirectedBreakpoint)new IdsvVariantContextBuilder(getContext(),
				new VariantContextBuilder()
					.chr("polyA")
					.start(1)
					.stop(1)
					.alleles("A", "A[polyA:10[")
					.make())
			.make()).getBreakendSummary();
		assertEquals(FWD, s1.direction);
		assertEquals(1, s1.start);
		assertEquals(1, s1.end);
		assertEquals(0, s1.referenceIndex);
		assertEquals(BWD, s1.direction2);
		assertEquals(10, s1.start2);
		assertEquals(10, s1.end2);
		assertEquals(0, s1.referenceIndex2);
	}
	@Test
	public void backward_breakend_should_match_vcf_call() {
		BreakpointSummary s1 = ((VariantContextDirectedBreakpoint)new IdsvVariantContextBuilder(getContext(),
				new VariantContextBuilder()
					.chr("polyA")
					.start(10)
					.stop(10)
					.alleles("A", "]polyA:1]A")
					.make())
			.make()).getBreakendSummary();
		assertEquals(BWD, s1.direction);
		assertEquals(10, s1.start);
		assertEquals(10, s1.end);
		assertEquals(0, s1.referenceIndex);
		assertEquals(FWD, s1.direction2);
		assertEquals(1, s1.start2);
		assertEquals(1, s1.end2);
		assertEquals(0, s1.referenceIndex2);
	}
	@Test
	public void breakpoint_width_should_be_determined_by_CIPOS() {
		BreakpointSummary s1 = ((VariantContextDirectedBreakpoint)new IdsvVariantContextBuilder(getContext(),
				new VariantContextBuilder()
					.chr("polyA")
					.start(2)
					.stop(2)
					.alleles("A", "A[polyA:10[")
					.attribute("CIPOS", new int[] { -1, 2 })
					.make())
			.make()).getBreakendSummary();
		assertEquals(1, s1.start);
		assertEquals(4, s1.end);
		assertEquals(10, s1.start2);
		assertEquals(10, s1.end2);
	}
	@Test
	public void remote_breakpoint_width_should_be_determined_by_CIRPOS() {
		BreakpointSummary s1 = ((VariantContextDirectedBreakpoint)new IdsvVariantContextBuilder(getContext(),
				new VariantContextBuilder()
					.chr("polyA")
					.start(2)
					.stop(2)
					.alleles("A", "A[polyA:10[")
					.attribute("CIRPOS", new int[] { -1, 2 })
					.make())
			.make()).getBreakendSummary();
		assertEquals(2, s1.start);
		assertEquals(2, s1.end);
		assertEquals(9, s1.start2);
		assertEquals(12, s1.end2);
	}
	@Test
	@Ignore
	public void percent_encoding_should_round_trip() {
		byte[] q = new byte[90];
		byte[] b = new byte[90];
		for (int i = 0; i < b.length; i++) {
			q[i] = (byte)i;
			b[i] = (byte)'A';
		}
		
		VariantContextDirectedEvidence e = (VariantContextDirectedEvidence) minimalBreakend()
			.breakend(new BreakendSummary(0, FWD, 1), b, q).make();
		e = (VariantContextDirectedEvidence)(IdsvVariantContext.create(getContext(), null, e));
		assertArrayEquals(q, e.getBreakendQuality());
	}
}
