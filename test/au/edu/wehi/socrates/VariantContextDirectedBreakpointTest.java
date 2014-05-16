package au.edu.wehi.socrates;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

import org.junit.Test;

public class VariantContextDirectedBreakpointTest extends TestHelper {
	@Test
	public void getBreakendSummary_should_handle_f_single_breakend() {
		BreakendSummary loc = new VariantContextDirectedBreakpoint(getContext(), minimalVariant().alleles("A", "A.").log10PError(-1).make()).getBreakendSummary();
		assertEquals(0, loc.referenceIndex);
		assertEquals(1, loc.start);
		assertEquals(1, loc.end);
		assertEquals(BreakendDirection.Forward, loc.direction);
	}
	@Test
	public void getBreakendSummary_should_handle_b_single_breakend() {
		BreakendSummary loc = new VariantContextDirectedBreakpoint(getContext(), minimalVariant().alleles("A", ".A").log10PError(-1).make()).getBreakendSummary();
		assertEquals(0, loc.referenceIndex);
		assertEquals(1, loc.start);
		assertEquals(1, loc.end);
		assertEquals(BreakendDirection.Backward, loc.direction);
	}
	@Test
	public void getEvidenceID_should_match_vcf_name() {
		VariantContextDirectedBreakpoint vc = new VariantContextDirectedBreakpoint(getContext(), minimalVariant().alleles("A", "AC[polyA:1[").id("test").make());
		assertEquals("test", vc.getEvidenceID());
	}
	@Test()
	public void should_allow_only_single_allele() {
		assertFalse(new VariantContextDirectedBreakpoint(getContext(), minimalVariant().alleles("A", "A.", "AA.").id("test").make()).isValid());
	}
	@Test
	public void getBreakpointSequence_should_match_f_single_alt_allele() {
		VariantContextDirectedBreakpoint vc = new VariantContextDirectedBreakpoint(getContext(), minimalVariant().alleles("A", "ACGT.").id("test").make());
		assertEquals("CGT", vc.getBreakpointSequenceString());
		assertArrayEquals(B(vc.getBreakpointSequenceString()), vc.getBreakpointSequence());
	}
	@Test
	public void getBreakpointSequence_should_match_b_single_alt_allele() {
		VariantContextDirectedBreakpoint vc = new VariantContextDirectedBreakpoint(getContext(), minimalVariant().alleles("A", ".CGTA").id("test").make());
		assertEquals("CGT", vc.getBreakpointSequenceString());
		assertArrayEquals(B(vc.getBreakpointSequenceString()), vc.getBreakpointSequence());
	}
	@Test
	public void getBreakpointSequence_should_return_untemplated_sequence_for_partnered_b_breakend() {
		VariantContextDirectedBreakpoint vc = new VariantContextDirectedBreakpoint(getContext(), minimalVariant().alleles("A", "[<polyA>[CGTA").id("test").make());
		assertEquals("CGT", vc.getBreakpointSequenceString());
		assertArrayEquals(B(vc.getBreakpointSequenceString()), vc.getBreakpointSequence());
	}
	@Test
	public void getBreakpointSequence_should_return_untemplated_sequence_for_partnered_f_breakend() {
		VariantContextDirectedBreakpoint vc = new VariantContextDirectedBreakpoint(getContext(), minimalVariant().alleles("A", "CGTA[<polyA>[").id("test").make());
		assertEquals("GTA", vc.getBreakpointSequenceString());
		assertArrayEquals(B(vc.getBreakpointSequenceString()), vc.getBreakpointSequence());
	}
	@Test
	public void getAnchorSequenceString_should_match_bases_corresponding_to_ref_positions() {
		VariantContextDirectedBreakpoint vc = new VariantContextDirectedBreakpoint(getContext(), minimalVariant().start(1).stop(4).alleles("ACGT", "AAAAT.").id("test").make());
		assertEquals("AAAA", vc.getAnchorSequenceString());
	}
	@Test
	public void getBreakendSummary_should_handle_partner_contig() {
		VariantContextDirectedBreakpoint vc = new VariantContextDirectedBreakpoint(getContext(), minimalVariant().start(1).stop(2).alleles("AA", "AAAAT[polyACGT:5[").id("test").make());
		BreakpointSummary loc = (BreakpointSummary)vc.getBreakendSummary();
		assertEquals(0, loc.referenceIndex);
		assertEquals(2, loc.start);
		assertEquals(2, loc.end);
		assertEquals(BreakendDirection.Forward, loc.direction);
		assertEquals(1, loc.referenceIndex2);
		assertEquals(5, loc.start2);
		assertEquals(5, loc.end2);
		assertEquals(BreakendDirection.Forward, loc.direction2);
	}
	@Test
	public void getBreakendSummary_should_handle_nonreference_partner_contig() {
		VariantContextDirectedBreakpoint vc = new VariantContextDirectedBreakpoint(getContext(), minimalVariant().start(1).stop(2).alleles("AA", "AAAAT[<contigNotInReference>:1[").id("test").make());
		BreakpointSummary loc = (BreakpointSummary)vc.getBreakendSummary();
		assertEquals(0, loc.referenceIndex);
		assertEquals(2, loc.start);
		assertEquals(2, loc.end);
		assertEquals(BreakendDirection.Forward, loc.direction);
		assertEquals(-1, loc.referenceIndex2);
		assertEquals(1, loc.start2);
		assertEquals(1, loc.end2);
		assertEquals(BreakendDirection.Forward, loc.direction2);
	}
	@Test
	public void getBreakendSummary_should_handle_ff_partner() {
		VariantContextDirectedBreakpoint vc = new VariantContextDirectedBreakpoint(getContext(), minimalVariant().start(1).stop(2).alleles("AA", "AAAAT[polyACGT:1234[").id("test").make());
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
		VariantContextDirectedBreakpoint vc = new VariantContextDirectedBreakpoint(getContext(), minimalVariant().start(1).stop(2).alleles("AA", "AAAAT]polyACGT:1234]").id("test").make());
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
		VariantContextDirectedBreakpoint vc = new VariantContextDirectedBreakpoint(getContext(), minimalVariant().start(1).stop(2).alleles("AA", "]polyACGT:1234]AAAAT").id("test").make());
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
		VariantContextDirectedBreakpoint vc = new VariantContextDirectedBreakpoint(getContext(), minimalVariant().start(1).stop(2).alleles("AA", "[polyACGT:1234[AAAAT").id("test").make());
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
}
