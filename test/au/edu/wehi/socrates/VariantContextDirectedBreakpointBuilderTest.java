package au.edu.wehi.socrates;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import org.junit.Test;

import au.edu.wehi.socrates.vcf.VcfAttributes;

public class VariantContextDirectedBreakpointBuilderTest extends TestHelper {
	@Test
	public void evidenceID_should_be_ID() {
		VariantContextDirectedBreakpointBuilder builder = new VariantContextDirectedBreakpointBuilder(getContext())
			.breakend(new BreakendSummary(0, FWD, 1, 1, null), null);
		builder.id("testID");
		VariantContextDirectedBreakpoint bp = builder.make();
		
		assertEquals("testID", bp.getID());
		assertEquals("testID", bp.getEvidenceID());
	}
	@Test
	public void null_reference_should_result_in_N_ref_call() {
		assertEquals("N", new VariantContextDirectedBreakpointBuilder(new ProcessingContext(null, getSequenceDictionary(), null))
			.breakend(new BreakendSummary(0, BreakendDirection.Forward, 1, 1, null), null).make().getReference().getDisplayString());
	}
	@Test
	public void should_lookup_reference_base() {
		VariantContextDirectedBreakpointBuilder builder = new VariantContextDirectedBreakpointBuilder(getContext());
		builder.breakend(new BreakendSummary(0, BreakendDirection.Forward, 1, 1, null), null);
		VariantContextDirectedBreakpoint bp = builder.make();
		
		assertEquals("A", bp.getReference().getDisplayString());
	}
	@Test
	public void should_round_trip_AssemblerProgram() {
		VariantContextDirectedBreakpoint dba = ReadEvidenceAssemblerUtil.create(getContext(), "test", 0, 1, BreakendDirection.Forward, B("A"), B("AA"), 1, 1, 0);
		assertEquals("test", dba.getAssemblerProgram());
		dba = new VariantContextDirectedBreakpoint(getContext(), new VariantContextBuilder(dba).make());
		assertEquals("test", dba.getAssemblerProgram());
	}
	@Test
	public void should_round_trip_getAssemblyConsensus() {
		VariantContextDirectedBreakpoint dba = ReadEvidenceAssemblerUtil.create(getContext(), "test", 0, 1, BreakendDirection.Forward, B("A"), B("AA"), 1, 1, 0);
		assertEquals("AA", dba.getAssemblyConsensus());
		dba = new VariantContextDirectedBreakpoint(getContext(), new VariantContextBuilder(dba).make());
		assertEquals("AA", dba.getAssemblyConsensus());
	}
	@Test
	public void getBreakpointSequence_should_get_untemplated_sequence() {
		VariantContextDirectedBreakpoint dba = ReadEvidenceAssemblerUtil.create(getContext(), "test", 0, 1, BreakendDirection.Forward, B("A"), B("AA"), 1, 1, 0);
		assertEquals("A", S(dba.getBreakpointSequence()));
		dba = new VariantContextDirectedBreakpoint(getContext(), new VariantContextBuilder(dba).make());
		assertEquals("A", S(dba.getBreakpointSequence()));
	}
	@Test
	public void should_match_variant_location_f() {
		VariantContextDirectedBreakpoint dba = ReadEvidenceAssemblerUtil.create(getContext(), "test", 0, 10, BreakendDirection.Forward, B("AT"), B("ATA"), 1, 1, 0);
		dba = new VariantContextDirectedBreakpoint(getContext(), new VariantContextBuilder(dba).make());
		assertEquals("polyA", dba.getChr());
		assertEquals(10, dba.getStart());
		assertEquals(10, dba.getEnd());
	}
	@Test
	public void should_match_variant_location_b() {
		VariantContextDirectedBreakpoint dba = ReadEvidenceAssemblerUtil.create(getContext(), "test", 0, 10, BreakendDirection.Backward, B("AT"), B("ATA"), 1, 1, 0);
		dba = new VariantContextDirectedBreakpoint(getContext(), new VariantContextBuilder(dba).make());
		assertEquals("polyA", dba.getChr());
		assertEquals(10, dba.getStart());
		assertEquals(10, dba.getEnd());
	}
	@Test
	public void id_should_be_based_on_assembler_position_direction() {
		VariantContextDirectedBreakpoint dba = ReadEvidenceAssemblerUtil.create(getContext(), "test", 0, 10, BreakendDirection.Backward, B("AT"), B("ATA"), 1, 1, 0);
		dba = new VariantContextDirectedBreakpoint(getContext(), new VariantContextBuilder(dba).make());
		assertEquals("test-polyA:10-b", dba.getID());
	}
	@Test
	public void phred_should_be_breakpoint_quality() {
		VariantContextDirectedBreakpoint dba = ReadEvidenceAssemblerUtil.create(getContext(), "test", 0, 10, BreakendDirection.Backward, B("AT"), B("ATA"), 1, 1, 7);
		dba = new VariantContextDirectedBreakpoint(getContext(), new VariantContextBuilder(dba).make());
		assertEquals(7, dba.getPhredScaledQual(), 0);
	}
	@Test
	public void should_generate_single_breakend_f() {
		VariantContextDirectedBreakpoint dba = ReadEvidenceAssemblerUtil.create(getContext(), "test", 0, 1, BreakendDirection.Forward, B("AT"), B("ATA"), 1, 1, 0);
		// ref base + breakpoint
		assertEquals("AAT.", dba.getAlternateAllele(0).getDisplayString());
		dba = new VariantContextDirectedBreakpoint(getContext(), new VariantContextBuilder(dba).make());
		assertEquals("AAT.", dba.getAlternateAllele(0).getDisplayString());
	}
	@Test
	public void should_generate_single_breakend_b() {
		// ref base + breakpoint
		VariantContextDirectedBreakpoint dba = ReadEvidenceAssemblerUtil.create(getContext(), "test", 0, 1, BreakendDirection.Backward, B("AT"), B("ATA"), 1, 1, 0);
		assertEquals(".ATA", dba.getAlternateAllele(0).getDisplayString());
		dba = new VariantContextDirectedBreakpoint(getContext(), new VariantContextBuilder(dba).make());
		assertEquals(".ATA", dba.getAlternateAllele(0).getDisplayString());
	}
	@Test
	public void should_set_realignment_failure_flag() {
		VariantContextDirectedBreakpointBuilder builder = new VariantContextDirectedBreakpointBuilder(getContext())
			.breakend(new BreakendSummary(0, FWD, 1, 1, null), null)
			.realignmentFailed();
		VariantContextDirectedBreakpoint bp = builder.make();
		assertTrue(bp.hasAttribute(VcfAttributes.REALIGNMENT_FAILURE.attribute()));
	}
	public VariantContextDirectedBreakpoint test_mated_breakend(BreakendDirection direction, boolean realignPositive, String bpString, String realignedCigar, String expectedAllele) {
		VariantContextDirectedBreakpoint dba = ReadEvidenceAssemblerUtil.create(getContext(), "test", 0, 1, direction, B(bpString), B(bpString), 1, 1, 0);
		dba = new VariantContextDirectedBreakpoint(getContext(), new VariantContextBuilder(dba).make());
		SAMRecord realigned = Read(1, 10, realignedCigar);
		realigned.setReadBases(realignPositive ? B(bpString) : B(SequenceUtil.reverseComplement(bpString)));
		realigned.setReadNegativeStrandFlag(!realignPositive);
		RealignedBreakpoint rbp = new RealignedBreakpoint(dba.getBreakendSummary(), realigned);
		VariantContextDirectedBreakpointBuilder builder = new VariantContextDirectedBreakpointBuilder(getContext(), dba)
			.breakend(rbp.getBreakpointSummary(), rbp.getInsertedSequence())
			.evidence(rbp.getBreakpointSummary().evidence);
		dba = builder.make();
		
		assertEquals(expectedAllele, dba.getAlternateAllele(0).getDisplayString());
		return dba;
	}
	@Test
	public void should_create_mated_breakend_if_realigned_ff() {
		// forward clip of CATCAT realigned to
		// 12345678901234567890
		//        .SMMMSS
		test_mated_breakend(BreakendDirection.Forward, true, "ATTTGC", "1S3M2S", "AA[polyACGT:10[");
	}
	@Test
	public void should_create_mated_breakend_if_realigned_fb() {
		// 12345678901234567890
		//         SMMMSS.
		test_mated_breakend(BreakendDirection.Forward, false, "ATTTGC", "1S3M2S", "AAT]polyACGT:12]");
	}
	@Test
	public void should_create_mated_breakend_if_realigned_bf() {
		// 12345678901234567890
		//         SMMMSS.
		test_mated_breakend(BreakendDirection.Backward, true, "ATTTGC", "1S3M2S", "]polyACGT:12]GCA");
	}
	@Test
	public void should_create_mated_breakend_if_realigned_bb() {
		// 12345678901234567890
		//        .SMMMSS
		test_mated_breakend(BreakendDirection.Backward, false, "ATTTGC", "1S3M2S", "[polyACGT:10[CA");
	}
	@Test
	public void breakend_should_be_vcf_sv_breakend() {
		VariantContextDirectedBreakpointBuilder builder = new VariantContextDirectedBreakpointBuilder(getContext())
			.breakend(new BreakendSummary(0, BreakendDirection.Forward, 1, 1, null), null);
		
		assertEquals("BND", builder.make().getAttributeAsString("SVTYPE", ""));
	}
	@Test
	public void breakpoint_should_be_vcf_sv_breakend() {
		VariantContextDirectedBreakpointBuilder builder = new VariantContextDirectedBreakpointBuilder(getContext())
			.breakpoint(new BreakpointSummary(0, BreakendDirection.Forward, 1, 1, 0, BreakendDirection.Forward, 1, 1, null), null);
		
		assertEquals("BND", builder.make().getAttributeAsString("SVTYPE", ""));
	}
	@Test
	public void assembly_should_set_assembly_properties() {
		VariantContextDirectedBreakpointBuilder builder = new VariantContextDirectedBreakpointBuilder(getContext())
			.breakend(new BreakendSummary(0, FWD, 1, 1, null), null)
			.assembly("programName", B("ACT"), new byte[] {1,2,3}, 7);
		VariantContextDirectedBreakpoint dba = builder.make();
		
		assertEquals("programName", dba.getAssemblerProgram());
		assertEquals("ACT", dba.getAssemblyConsensus());
		assertEquals(7, dba.getAssemblyQuality(), 0);
	}
	@Test
	public void evidence_should_round_trip() {
		for (VcfAttributes a : VcfAttributes.evidenceValues()) {
			EvidenceMetrics em = new EvidenceMetrics();
			em.set(a, 1);
			VariantContextDirectedBreakpointBuilder builder = new VariantContextDirectedBreakpointBuilder(getContext())
				.breakend(new BreakendSummary(0, FWD, 1, 1, null), null)
				.evidence(em);
			builder = new VariantContextDirectedBreakpointBuilder(getContext(), builder.make());
			
			assertEquals(1, builder.make().getBreakendSummary().evidence.get(a));
		}
	}
	@Test
	public void evidence_score_should_match_vcf_phred_score() {
		EvidenceMetrics m = new EvidenceMetrics() {
			public double getScore() { return 10; };
		};
		VariantContextDirectedBreakpointBuilder builder = new VariantContextDirectedBreakpointBuilder(getContext())
			.breakend(new BreakendSummary(0, FWD, 1, 1, null), null)
			.evidence(m);
		assertEquals(10, builder.make().getPhredScaledQual(), 0);
	}
	@Test
	public void should_round_trip_inexact_breakend() {
		VariantContextDirectedBreakpointBuilder builder = new VariantContextDirectedBreakpointBuilder(getContext());
		builder.breakend(new BreakendSummary(1, FWD, 2, 4, null), null);
		VariantContextDirectedBreakpoint v = new VariantContextDirectedBreakpointBuilder(getContext(), builder.make()).make();
		assertEquals(1, v.getBreakendSummary().referenceIndex);
		assertEquals(FWD, v.getBreakendSummary().direction);
		assertEquals(2, v.getBreakendSummary().start);
		assertEquals(4, v.getBreakendSummary().end);
	}
	@Test
	public void should_write_breakpoint_in_vcf41_mode() {
		ProcessingContext context = getContext();
		context.setVcf41Mode(true);
		VariantContextDirectedBreakpointBuilder builder = new VariantContextDirectedBreakpointBuilder(context);
		builder.breakend(new BreakendSummary(0,  FWD,  1,  1, null), "ACGT");
		VariantContextDirectedBreakpoint vc = builder.make();
		assertEquals(-1, vc.getAlternateAlleles().get(0).getDisplayString().indexOf("."));
	}
	@Test
	public void vcf41_breakend_should_round_trip() {
		ProcessingContext context = getContext();
		context.setVcf41Mode(true);
		VariantContextDirectedBreakpointBuilder builder = new VariantContextDirectedBreakpointBuilder(context);
		builder.breakend(new BreakendSummary(0,  FWD,  1,  2, null), "ACGT");
		VariantContextDirectedBreakpoint v = new VariantContextDirectedBreakpointBuilder(context, builder.make()).make();
		assertEquals(0, v.getBreakendSummary().referenceIndex);
		assertEquals(FWD, v.getBreakendSummary().direction);
		assertEquals(1, v.getBreakendSummary().start);
		assertEquals(2, v.getBreakendSummary().end);
		
		builder = new VariantContextDirectedBreakpointBuilder(context);
		builder.breakend(new BreakendSummary(0,  BWD,  1,  2, null), "ACGT");
		v = new VariantContextDirectedBreakpointBuilder(context, builder.make()).make();
		assertEquals(0, v.getBreakendSummary().referenceIndex);
		assertEquals(BWD, v.getBreakendSummary().direction);
		assertEquals(1, v.getBreakendSummary().start);
		assertEquals(2, v.getBreakendSummary().end);
	}
}