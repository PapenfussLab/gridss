package au.edu.wehi.socrates;

import net.sf.samtools.SAMRecord;

import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.junit.Test;
import org.junit.rules.ExpectedException;

import static org.junit.Assert.*;

public class DirectedBreakpointAssemblyTest extends TestHelper {
	@Test
	public void null_reference_should_result_in_N_ref_call() {
		assertEquals("N", DirectedBreakpointAssembly.create(getSequenceDictionary(), null, "test", 0, 1, BreakpointDirection.Forward, B("A"), B("AA"), 1, 0).getReference().getBaseString());
	}
	@Test
	public void should_lookup_reference_base() {
		assertEquals("A", DirectedBreakpointAssembly.create(getSequenceDictionary(), SMALL_FA, "test", 0, 1, BreakpointDirection.Forward, B("A"), B("AA"), 1, 0).getReference().getBaseString());
	}
	@Test
	public void should_round_trip_AssemblerProgram() {
		DirectedBreakpointAssembly dba = DirectedBreakpointAssembly.create(getSequenceDictionary(), SMALL_FA, "test", 0, 1, BreakpointDirection.Forward, B("A"), B("AA"), 1, 0);
		assertEquals("test", dba.getAssemblerProgram());
		dba = new DirectedBreakpointAssembly(getSequenceDictionary(), new VariantContextBuilder(dba).make());
		assertEquals("test", dba.getAssemblerProgram());
	}
	@Test
	public void should_round_trip_getAssemblyConsensus() {
		DirectedBreakpointAssembly dba = DirectedBreakpointAssembly.create(getSequenceDictionary(), SMALL_FA, "test", 0, 1, BreakpointDirection.Forward, B("A"), B("AA"), 1, 0);
		assertEquals("test", dba.getAssemblyConsensus());
		dba = new DirectedBreakpointAssembly(getSequenceDictionary(), new VariantContextBuilder(dba).make());
		assertEquals("test", dba.getAssemblyConsensus());
	}
	@Test
	public void should_round_trip_getAssemblyReadCount() {
		DirectedBreakpointAssembly dba = DirectedBreakpointAssembly.create(getSequenceDictionary(), SMALL_FA, "test", 0, 1, BreakpointDirection.Forward, B("A"), B("AA"), 7, 0);
		assertEquals(7, dba.getConsensusReadCount());
		dba = new DirectedBreakpointAssembly(getSequenceDictionary(), new VariantContextBuilder(dba).make());
		assertEquals(7, dba.getConsensusReadCount());
	}
	@Test
	public void should_match_variant_location_f() {
		DirectedBreakpointAssembly dba = DirectedBreakpointAssembly.create(getSequenceDictionary(), SMALL_FA, "test", 0, 10, BreakpointDirection.Forward, B("AT"), B("ATA"), 1, 0);
		dba = new DirectedBreakpointAssembly(getSequenceDictionary(), new VariantContextBuilder(dba).make());
		assertEquals("polyA", dba.getChr());
		assertEquals(10, dba.getStart());
		assertEquals(10, dba.getEnd());
	}
	@Test
	public void should_match_variant_location_b() {
		DirectedBreakpointAssembly dba = DirectedBreakpointAssembly.create(getSequenceDictionary(), SMALL_FA, "test", 0, 10, BreakpointDirection.Backward, B("AT"), B("ATA"), 1, 0);
		dba = new DirectedBreakpointAssembly(getSequenceDictionary(), new VariantContextBuilder(dba).make());
		assertEquals("polyA", dba.getChr());
		assertEquals(10, dba.getStart());
		assertEquals(10, dba.getEnd());
	}
	@Test
	public void id_should_be_based_on_assembler_position_direction() {
		DirectedBreakpointAssembly dba = DirectedBreakpointAssembly.create(getSequenceDictionary(), SMALL_FA, "test", 0, 10, BreakpointDirection.Backward, B("AT"), B("ATA"), 1, 0);
		dba = new DirectedBreakpointAssembly(getSequenceDictionary(), new VariantContextBuilder(dba).make());
		assertEquals("test-polyA:10-b", dba.getID());
	}
	@Test
	public void phred_should_be_breakpoint_quality() {
		DirectedBreakpointAssembly dba = DirectedBreakpointAssembly.create(getSequenceDictionary(), SMALL_FA, "test", 0, 10, BreakpointDirection.Backward, B("AT"), B("ATA"), 1, 7);
		dba = new DirectedBreakpointAssembly(getSequenceDictionary(), new VariantContextBuilder(dba).make());
		assertEquals(7, dba.getLog10PError(), 0);
	}
	@Test
	public void should_generate_single_breakend_f() {
		DirectedBreakpointAssembly dba = DirectedBreakpointAssembly.create(getSequenceDictionary(), SMALL_FA, "test", 0, 1, BreakpointDirection.Forward, B("AT"), B("ATA"), 1, 0);
		// ref base + breakpoint
		assertEquals("AAT.", dba.getAlternateAllele(0).getBaseString());
		dba = new DirectedBreakpointAssembly(getSequenceDictionary(), new VariantContextBuilder(dba).make());
		assertEquals("AAT.", dba.getAlternateAllele(0).getBaseString());
	}
	@Test
	public void should_generate_single_breakend_b() {
		// ref base + breakpoint
		DirectedBreakpointAssembly dba = DirectedBreakpointAssembly.create(getSequenceDictionary(), SMALL_FA, "test", 0, 1, BreakpointDirection.Backward, B("AT"), B("ATA"), 1, 0);
		assertEquals(".ATA", dba.getAlternateAllele(0).getBaseString());
		dba = new DirectedBreakpointAssembly(getSequenceDictionary(), new VariantContextBuilder(dba).make());
		assertEquals(".ATA", dba.getAlternateAllele(0).getBaseString());
	}
	@Test
	public void should_not_create_mated_breakend_if_realigned_unmapped() {
		DirectedBreakpointAssembly dba = DirectedBreakpointAssembly.create(getSequenceDictionary(), SMALL_FA, "test", 0, 1, BreakpointDirection.Backward, B("AT"), B("ATA"), 1, 0);
		dba = new DirectedBreakpointAssembly(getSequenceDictionary(), new VariantContextBuilder(dba).make());
		dba = DirectedBreakpointAssembly.create(dba, Unmapped(2));
		assertEquals(".ATA", dba.getAlternateAllele(0).getBaseString());
	}
	public DirectedBreakpointAssembly test_mated_breakend(BreakpointDirection direction, boolean mateNegative, String bpString, String realignedCigar, String expectedAllele) {
		DirectedBreakpointAssembly dba = DirectedBreakpointAssembly.create(getSequenceDictionary(), null, "test", 0, 1, direction, B(bpString), B(bpString), 1, 0);
		dba = new DirectedBreakpointAssembly(getSequenceDictionary(), new VariantContextBuilder(dba).make());
		SAMRecord realigned = Read(1, 10, realignedCigar);
		realigned.setReadNegativeStrandFlag(mateNegative);
		dba = DirectedBreakpointAssembly.create(dba, realigned);
		assertEquals(expectedAllele, dba.getAlternateAllele(0).getBaseString());
		return dba;
	}
	@Test
	public void should_create_mated_breakend_if_realigned_ff() {
		test_mated_breakend(BreakpointDirection.Forward, true, "CATCAT", "1S3M2S", "NC[polyACGT:10[");
	}
	@Test
	public void should_create_mated_breakend_if_realigned_fb() {
		test_mated_breakend(BreakpointDirection.Forward, false, "CATCAT", "1S3M2S", "NCA]polyACGT:10]");
	}
	@Test
	public void should_create_mated_breakend_if_realigned_bf() {
		test_mated_breakend(BreakpointDirection.Backward, true, "CATCAT", "1S3M2S", "]polyACGT:10]ATN");
	}
	@Test
	public void should_create_mated_breakend_if_realigned_bb() {
		test_mated_breakend(BreakpointDirection.Backward, false, "CATCAT", "1S3M2S", "[polyACGT:10[TN");
	}
}
