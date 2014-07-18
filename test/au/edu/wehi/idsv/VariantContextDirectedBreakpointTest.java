package au.edu.wehi.idsv;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import org.junit.Ignore;
import org.junit.Test;

public class VariantContextDirectedBreakpointTest extends TestHelper {
	@Test
	public void getBreakendSummary_should_handle_f_single_breakend() {
		BreakendSummary loc = new VariantContextDirectedBreakpoint(getContext(), AES(), minimalVariant().alleles("A", "A.").log10PError(-1).make()).getBreakendSummary();
		assertEquals(0, loc.referenceIndex);
		assertEquals(1, loc.start);
		assertEquals(1, loc.end);
		assertEquals(BreakendDirection.Forward, loc.direction);
	}
	@Test
	public void getBreakendSummary_should_handle_b_single_breakend() {
		BreakendSummary loc = new VariantContextDirectedBreakpoint(getContext(), AES(), minimalVariant().alleles("A", ".A").log10PError(-1).make()).getBreakendSummary();
		assertEquals(0, loc.referenceIndex);
		assertEquals(1, loc.start);
		assertEquals(1, loc.end);
		assertEquals(BreakendDirection.Backward, loc.direction);
	}
	@Test
	public void getEvidenceID_should_match_vcf_name() {
		VariantContextDirectedBreakpoint vc = new VariantContextDirectedBreakpoint(getContext(), AES(), minimalVariant().alleles("A", "AC[polyA:1[").id("test").make());
		assertEquals("test", vc.getEvidenceID());
	}
	@Test()
	public void should_allow_only_single_allele() {
		assertFalse(new VariantContextDirectedBreakpoint(getContext(), AES(), minimalVariant().alleles("A", "A.", "AA.").id("test").make()).isValid());
	}
	@Test
	public void getBreakpointSequence_should_match_f_single_alt_allele() {
		VariantContextDirectedBreakpoint vc = new VariantContextDirectedBreakpoint(getContext(), AES(), minimalVariant().alleles("A", "ACGT.").id("test").make());
		assertEquals("CGT", vc.getBreakpointSequenceString());
		assertArrayEquals(B(vc.getBreakpointSequenceString()), vc.getBreakendSequence());
	}
	@Test
	public void getBreakpointSequence_should_match_b_single_alt_allele() {
		VariantContextDirectedBreakpoint vc = new VariantContextDirectedBreakpoint(getContext(), AES(), minimalVariant().alleles("A", ".CGTA").id("test").make());
		assertEquals("CGT", vc.getBreakpointSequenceString());
		assertArrayEquals(B(vc.getBreakpointSequenceString()), vc.getBreakendSequence());
	}
	@Test
	public void getBreakpointSequence_should_return_untemplated_sequence_for_partnered_b_breakend() {
		VariantContextDirectedBreakpoint vc = new VariantContextDirectedBreakpoint(getContext(), AES(), minimalVariant().alleles("A", "[polyA[CGTA").id("test").make());
		assertEquals("CGT", vc.getBreakpointSequenceString());
		assertArrayEquals(B(vc.getBreakpointSequenceString()), vc.getBreakendSequence());
	}
	@Test
	public void getBreakpointSequence_should_return_untemplated_sequence_for_partnered_f_breakend() {
		VariantContextDirectedBreakpoint vc = new VariantContextDirectedBreakpoint(getContext(), AES(), minimalVariant().alleles("A", "CGTA[polyA[").id("test").make());
		assertEquals("GTA", vc.getBreakpointSequenceString());
		assertArrayEquals(B(vc.getBreakpointSequenceString()), vc.getBreakendSequence());
	}
	@Test
	public void getAnchorSequenceString_should_match_bases_corresponding_to_ref_positions() {
		VariantContextDirectedBreakpoint vc = new VariantContextDirectedBreakpoint(getContext(), AES(), minimalVariant().start(1).stop(4).alleles("ACGT", "AAAAT.").id("test").make());
		assertEquals("AAAA", vc.getAnchorSequenceString());
	}
	@Test
	public void getBreakendSummary_should_handle_partner_contig() {
		VariantContextDirectedBreakpoint vc = new VariantContextDirectedBreakpoint(getContext(), AES(), minimalVariant().start(1).stop(2).alleles("AA", "AAAAT[polyACGT:5[").id("test").make());
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
		ProcessingContext pc = getContext();
		int seqCount = pc.getDictionary().size();
		new VariantContextDirectedBreakpoint(pc, AES(), minimalVariant().chr("<contigNotInReference>").make());
		assertEquals(seqCount + 1, pc.getDictionary().size());
	}
	@Test
	public void should_extend_process_context_sequence_dictionary_when_encountering_named_contig_partner() {
		ProcessingContext pc = getContext();
		int seqCount = pc.getDictionary().size();
		new VariantContextDirectedBreakpoint(pc, AES(), minimalVariant().start(1).stop(2).alleles("AA", "AAAAT[<contigNotInReference>:1[").id("test").make());
		assertEquals(seqCount + 1, pc.getDictionary().size());
	}
	@Test
	public void should_not_extend_process_context_sequence_dictionary_when_encountering_vcf41_placeholder() {
		ProcessingContext pc = getContext();
		int seqCount = pc.getDictionary().size();
		new VariantContextDirectedBreakpoint(pc, AES(), minimalVariant().start(1).stop(2).alleles("AA", "AAAAT[<IDSV_PLACEHOLDER>[").id("test").make());
		assertEquals(seqCount, pc.getDictionary().size());
	}
	@Test
	public void getBreakendSummary_should_handle_ff_partner() {
		VariantContextDirectedBreakpoint vc = new VariantContextDirectedBreakpoint(getContext(), AES(), minimalVariant().start(1).stop(2).alleles("AA", "AAAAT]polyACGT:1234]").id("test").make());
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
		VariantContextDirectedBreakpoint vc = new VariantContextDirectedBreakpoint(getContext(), AES(), minimalVariant().start(1).stop(2).alleles("AA", "AAAAT[polyACGT:1234[").id("test").make());
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
		VariantContextDirectedBreakpoint vc = new VariantContextDirectedBreakpoint(getContext(), AES(), minimalVariant().start(1).stop(2).alleles("AA", "[polyACGT:1234[AAAAT").id("test").make());
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
		VariantContextDirectedBreakpoint vc = new VariantContextDirectedBreakpoint(getContext(), AES(), minimalVariant().start(1).stop(2).alleles("AA", "]polyACGT:1234]AAAAT").id("test").make());
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
	@Test
	public void getBreakendSummary_should_parse_vcf41_compatability_breakpoint_as_breakend() {
		// even in 4.2 mode, we should happily parse our backward compatible serialisation as a breakend
		VariantContextDirectedBreakpoint vc = new VariantContextDirectedBreakpoint(getContext(), AES(), minimalVariant().start(1).stop(1).alleles("A", "A[<IDSV_PLACEHOLDER>[").id("test").make());
		assertEquals(vc.getBreakendSummary().getClass(), BreakendSummary.class);
	}
	@Test
	public void matching_breakpoints_should_call_same_evidence() {
		BreakpointSummary s1 = (BreakpointSummary)new VariantContextDirectedBreakpointBuilder(getContext(), AES(),
				new VariantContextBuilder()
					.chr("polyA")
					.start(1)
					.stop(1)
					.alleles("A", "A[polyA:10[")
					.make())
			.make().getBreakendSummary();
		BreakpointSummary s2 = (BreakpointSummary)new VariantContextDirectedBreakpointBuilder(getContext(), AES(),
				new VariantContextBuilder()
					.chr("polyA")
					.start(10)
					.stop(10)
					.alleles("A", "]polyA:1]A")
					.make())
			.make().getBreakendSummary();
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
		BreakpointSummary s1 = (BreakpointSummary)new VariantContextDirectedBreakpointBuilder(getContext(), AES(),
				new VariantContextBuilder()
					.chr("polyA")
					.start(1)
					.stop(1)
					.alleles("A", "A[polyA:10[")
					.make())
			.make().getBreakendSummary();
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
		BreakpointSummary s1 = (BreakpointSummary)new VariantContextDirectedBreakpointBuilder(getContext(), AES(),
				new VariantContextBuilder()
					.chr("polyA")
					.start(10)
					.stop(10)
					.alleles("A", "]polyA:1]A")
					.make())
			.make().getBreakendSummary();
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
		BreakpointSummary s1 = (BreakpointSummary)new VariantContextDirectedBreakpointBuilder(getContext(), AES(),
				new VariantContextBuilder()
					.chr("polyA")
					.start(2)
					.stop(2)
					.alleles("A", "A[polyA:10[")
					.attribute("CIPOS", new int[] { -1, 2 })
					.make())
			.make().getBreakendSummary();
		assertEquals(1, s1.start);
		assertEquals(4, s1.end);
		assertEquals(10, s1.start2);
		assertEquals(10, s1.end2);
	}
	@Test
	public void remote_breakpoint_width_should_be_determined_by_CIRPOS() {
		BreakpointSummary s1 = (BreakpointSummary)new VariantContextDirectedBreakpointBuilder(getContext(), AES(),
				new VariantContextBuilder()
					.chr("polyA")
					.start(2)
					.stop(2)
					.alleles("A", "A[polyA:10[")
					.attribute("CIRPOS", new int[] { -1, 2 })
					.make())
			.make().getBreakendSummary();
		assertEquals(2, s1.start);
		assertEquals(2, s1.end);
		assertEquals(9, s1.start2);
		assertEquals(12, s1.end2);
	}
}
