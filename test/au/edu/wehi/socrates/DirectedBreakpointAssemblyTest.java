package au.edu.wehi.socrates;

import static org.junit.Assert.*;
import htsjdk.samtools.SAMRecord;

import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.junit.Test;

import au.edu.wehi.socrates.vcf.VcfSvConstants;

public class DirectedBreakpointAssemblyTest extends TestHelper {
	@Test
	public void null_reference_should_result_in_N_ref_call() {
		assertEquals("N", DirectedBreakpointAssembly.create(new ProcessingContext(null, getSequenceDictionary(), null), "test", 0, 1, BreakendDirection.Forward, B("A"), B("AA"), 1, 0).getReference().getDisplayString());
	}
	@Test
	public void should_lookup_reference_base() {
		assertEquals("A", DirectedBreakpointAssembly.create(getContext(), "test", 0, 1, BreakendDirection.Forward, B("A"), B("AA"), 1, 0).getReference().getDisplayString());
	}
	@Test
	public void anchor_should_use_reference_base_not_assembly_base() {
		String alt = DirectedBreakpointAssembly.create(getContext(), "test", 0, 1, BreakendDirection.Forward, B("A"), B("CC"), 1, 0).getAlternateAllele(0).getDisplayString();
		assertEquals('A', alt.charAt(0));
	}
	@Test
	public void should_round_trip_AssemblerProgram() {
		DirectedBreakpointAssembly dba = DirectedBreakpointAssembly.create(getContext(), "test", 0, 1, BreakendDirection.Forward, B("A"), B("AA"), 1, 0);
		assertEquals("test", dba.getAssemblerProgram());
		dba = new DirectedBreakpointAssembly(getContext(), new VariantContextBuilder(dba).make());
		assertEquals("test", dba.getAssemblerProgram());
	}
	@Test
	public void getBreakpointSequence_should_get_untemplate_sequence() {
		DirectedBreakpointAssembly dba = DirectedBreakpointAssembly.create(getContext(), "test", 0, 1, BreakendDirection.Forward, B("A"), B("AA"), 1, 0);
		assertEquals("A", S(dba.getBreakpointSequence()));
		dba = new DirectedBreakpointAssembly(getContext(), new VariantContextBuilder(dba).make());
		assertEquals("A", S(dba.getBreakpointSequence()));
	}
	@Test
	public void should_round_trip_getAssemblyConsensus() {
		DirectedBreakpointAssembly dba = DirectedBreakpointAssembly.create(getContext(), "test", 0, 1, BreakendDirection.Forward, B("A"), B("AA"), 1, 0);
		assertEquals("AA", dba.getAssemblyConsensus());
		dba = new DirectedBreakpointAssembly(getContext(), new VariantContextBuilder(dba).make());
		assertEquals("AA", dba.getAssemblyConsensus());
	}
	@Test
	public void should_round_trip_getAssemblyReadCount() {
		DirectedBreakpointAssembly dba = DirectedBreakpointAssembly.create(getContext(), "test", 0, 1, BreakendDirection.Forward, B("A"), B("AA"), 7, 0);
		assertEquals(7, dba.getConsensusReadCount());
		dba = new DirectedBreakpointAssembly(getContext(), new VariantContextBuilder(dba).make());
		assertEquals(7, dba.getConsensusReadCount());
	}
	@Test
	public void should_match_variant_location_f() {
		DirectedBreakpointAssembly dba = DirectedBreakpointAssembly.create(getContext(), "test", 0, 10, BreakendDirection.Forward, B("AT"), B("ATA"), 1, 0);
		dba = new DirectedBreakpointAssembly(getContext(), new VariantContextBuilder(dba).make());
		assertEquals("polyA", dba.getChr());
		assertEquals(10, dba.getStart());
		assertEquals(10, dba.getEnd());
	}
	@Test
	public void should_match_variant_location_b() {
		DirectedBreakpointAssembly dba = DirectedBreakpointAssembly.create(getContext(), "test", 0, 10, BreakendDirection.Backward, B("AT"), B("ATA"), 1, 0);
		dba = new DirectedBreakpointAssembly(getContext(), new VariantContextBuilder(dba).make());
		assertEquals("polyA", dba.getChr());
		assertEquals(10, dba.getStart());
		assertEquals(10, dba.getEnd());
	}
	@Test
	public void id_should_be_based_on_assembler_position_direction() {
		DirectedBreakpointAssembly dba = DirectedBreakpointAssembly.create(getContext(), "test", 0, 10, BreakendDirection.Backward, B("AT"), B("ATA"), 1, 0);
		dba = new DirectedBreakpointAssembly(getContext(), new VariantContextBuilder(dba).make());
		assertEquals("test-polyA:10-b", dba.getID());
	}
	@Test
	public void phred_should_be_breakpoint_quality() {
		DirectedBreakpointAssembly dba = DirectedBreakpointAssembly.create(getContext(), "test", 0, 10, BreakendDirection.Backward, B("AT"), B("ATA"), 1, 7);
		dba = new DirectedBreakpointAssembly(getContext(), new VariantContextBuilder(dba).make());
		assertEquals(7, dba.getPhredScaledQual(), 0);
	}
	@Test
	public void should_generate_single_breakend_f() {
		DirectedBreakpointAssembly dba = DirectedBreakpointAssembly.create(getContext(), "test", 0, 1, BreakendDirection.Forward, B("AT"), B("ATA"), 1, 0);
		// ref base + breakpoint
		assertEquals("AAT.", dba.getAlternateAllele(0).getDisplayString());
		dba = new DirectedBreakpointAssembly(getContext(), new VariantContextBuilder(dba).make());
		assertEquals("AAT.", dba.getAlternateAllele(0).getDisplayString());
	}
	@Test
	public void should_generate_single_breakend_b() {
		// ref base + breakpoint
		DirectedBreakpointAssembly dba = DirectedBreakpointAssembly.create(getContext(), "test", 0, 1, BreakendDirection.Backward, B("AT"), B("ATA"), 1, 0);
		assertEquals(".ATA", dba.getAlternateAllele(0).getDisplayString());
		dba = new DirectedBreakpointAssembly(getContext(), new VariantContextBuilder(dba).make());
		assertEquals(".ATA", dba.getAlternateAllele(0).getDisplayString());
	}
	@Test
	public void should_not_create_mated_breakend_if_realigned_unmapped() {
		DirectedBreakpointAssembly dba = DirectedBreakpointAssembly.create(getContext(), "test", 0, 1, BreakendDirection.Backward, B("AT"), B("ATA"), 1, 0);
		dba = new DirectedBreakpointAssembly(getContext(), new VariantContextBuilder(dba).make());
		dba = DirectedBreakpointAssembly.create(dba, Unmapped(2));
		assertEquals(".ATA", dba.getAlternateAllele(0).getDisplayString());
	}
	public DirectedBreakpointAssembly test_mated_breakend(BreakendDirection direction, boolean realignPositive, String bpString, String realignedCigar, String expectedAllele) {
		DirectedBreakpointAssembly dba = DirectedBreakpointAssembly.create(getContext(), "test", 0, 1, direction, B(bpString), B(bpString), 1, 0);
		dba = new DirectedBreakpointAssembly(getContext(), new VariantContextBuilder(dba).make());
		SAMRecord realigned = Read(1, 10, realignedCigar);
		realigned.setReadNegativeStrandFlag(!realignPositive);
		dba = DirectedBreakpointAssembly.create(dba, realigned);
		assertEquals(expectedAllele, dba.getAlternateAllele(0).getDisplayString());
		return dba;
	}
	@Test
	public void should_create_mated_breakend_if_realigned_ff() {
		// forward clip of CATCAT realigned to
		// 12345678901234567890
		//        .SMMMSS
		//         CATCAT
		test_mated_breakend(BreakendDirection.Forward, true, "CGTCGT", "1S3M2S", "AC[polyACGT:10[");
	}
	@Test
	public void should_create_mated_breakend_if_realigned_fb() {
		// 12345678901234567890
		//         SMMMSS.
		//         ATGATG
		test_mated_breakend(BreakendDirection.Forward, false, "CGTCGT", "1S3M2S", "ACG]polyACGT:12]");
	}
	@Test
	public void should_create_mated_breakend_if_realigned_bf() {
		// 12345678901234567890
		//         SMMMSS.
		//         CATCAT
		test_mated_breakend(BreakendDirection.Backward, true, "CGTCGT", "1S3M2S", "]polyACGT:12]GTA");
	}
	@Test
	public void should_create_mated_breakend_if_realigned_bb() {
		// 12345678901234567890
		//        .SMMMSS
		//         ATGATG
		test_mated_breakend(BreakendDirection.Backward, false, "CGTCGT", "1S3M2S", "[polyACGT:10[TA");
	}
	@Test
	public void should_be_vcf_sv_breakend() {
		assertEquals("BND", DirectedBreakpointAssembly.create(getContext(), "test", 0, 1, BreakendDirection.Forward, B("A"), B("AA"), 1, 0)
			.getAttributeAsString(VcfSvConstants.SV_TYPE_KEY, null));
		assertEquals("BND", DirectedBreakpointAssembly.create(DirectedBreakpointAssembly.create(getContext(), "test", 0, 1, BreakendDirection.Forward, B("A"), B("AA"), 1, 0), Read(0, 1, "10M"))
				.getAttributeAsString(VcfSvConstants.SV_TYPE_KEY, null));
	}
	@Test
	public void should_assembly_base_qualities() {
		DirectedBreakpointAssembly dba = DirectedBreakpointAssembly.create(getContext(), "test", 0, 1, BreakendDirection.Forward, B("TTT"), new byte[] { 1,2,3 }, B("ATTT"), new byte[] { 1,2,3,4}, 1, 0);
		assertArrayEquals(new byte[] {1, 2, 3}, dba.getBreakpointQuality());
	}
}
