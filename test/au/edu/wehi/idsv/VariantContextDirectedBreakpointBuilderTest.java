package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.Header;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import org.junit.Test;

import au.edu.wehi.idsv.vcf.VcfAttributes;

public class VariantContextDirectedBreakpointBuilderTest extends TestHelper {
	@Test
	public void evidenceID_should_be_ID() {
		VariantContextDirectedBreakpointBuilder builder = new VariantContextDirectedBreakpointBuilder(getContext(), AES())
			.breakend(new BreakendSummary(0, FWD, 1, 1, null), null);
		builder.id("testID");
		VariantContextDirectedBreakpoint bp = builder.make();
		
		assertEquals("testID", bp.getID());
		assertEquals("testID", bp.getEvidenceID());
	}
	@Test
	public void should_lookup_reference_base() {
		VariantContextDirectedBreakpointBuilder builder = new VariantContextDirectedBreakpointBuilder(getContext(), AES());
		builder.breakend(new BreakendSummary(0, BreakendDirection.Forward, 1, 1, null), null);
		VariantContextDirectedBreakpoint bp = builder.make();
		
		assertEquals("A", bp.getReference().getDisplayString());
	}
	@Test
	public void should_round_trip_AssemblerProgram() {
		VariantContextDirectedBreakpoint dba = AE();
		assertEquals("testAssembler", dba.getAssemblerProgram());
		dba = new VariantContextDirectedBreakpoint(getContext(), AES(), new VariantContextBuilder(dba).make());
		assertEquals("testAssembler", dba.getAssemblerProgram());
	}
	@Test
	public void should_round_trip_getAssemblyConsensus() {
		VariantContextDirectedBreakpoint dba = AE();
		assertEquals("ATT", dba.getAssemblyConsensus());
		dba = new VariantContextDirectedBreakpoint(getContext(), AES(), new VariantContextBuilder(dba).make());
		assertEquals("ATT", dba.getAssemblyConsensus());
	}
	@Test
	public void getBreakpointSequence_should_get_untemplated_sequence() {
		VariantContextDirectedBreakpoint dba = AE();
		assertEquals("TT", S(dba.getBreakpointSequence()));
		dba = new VariantContextDirectedBreakpoint(getContext(), AES(), new VariantContextBuilder(dba).make());
		assertEquals("TT", S(dba.getBreakpointSequence()));
	}
	@Test
	public void should_match_variant_location_f() {
		VariantContextDirectedBreakpoint dba = AB().referenceAnchor(0, 10).makeVariant(); 
		dba = new VariantContextDirectedBreakpoint(getContext(), AES(), new VariantContextBuilder(dba).make());
		assertEquals("polyA", dba.getChr());
		assertEquals(10, dba.getStart());
		assertEquals(10, dba.getEnd());
	}
	@Test
	public void should_match_variant_location_b() {
		VariantContextDirectedBreakpoint dba = AB().referenceAnchor(0, 10).direction(BWD).makeVariant();
		dba = new VariantContextDirectedBreakpoint(getContext(), AES(), new VariantContextBuilder(dba).make());
		assertEquals("polyA", dba.getChr());
		assertEquals(10, dba.getStart());
		assertEquals(10, dba.getEnd());
	}
	@Test
	public void id_should_be_based_on_assembler_position_direction() {
		VariantContextDirectedBreakpoint dba = AE();
		dba = new VariantContextDirectedBreakpoint(getContext(), AES(), new VariantContextBuilder(dba).make());
		assertTrue(dba.getID().startsWith("testAssembler-polyA:1-f"));
	}
	@Test
	public void phred_should_be_evidence_score() {
		VariantContextDirectedBreakpoint dba = AE();
		dba = new VariantContextDirectedBreakpoint(getContext(), AES(), new VariantContextBuilder(dba).make());
		assertEquals(dba.getBreakendSummary().evidence.getScore(), dba.getPhredScaledQual(), 0);
	}
	@Test
	public void should_generate_single_breakend_f() {
		VariantContextDirectedBreakpoint dba = AB().assemblyBases(B("NNGT")).anchorLength(2).makeVariant(); 
		// ref base + breakpoint
		assertEquals("AGT.", dba.getAlternateAllele(0).getDisplayString());
		dba = new VariantContextDirectedBreakpoint(getContext(), AES(), new VariantContextBuilder(dba).make());
		assertEquals("AGT.", dba.getAlternateAllele(0).getDisplayString());
	}
	@Test
	public void should_generate_single_breakend_b() {
		// ref base + breakpoint
		VariantContextDirectedBreakpoint dba = AB().assemblyBases(B("GTNN")).anchorLength(2).direction(BWD).makeVariant();
		assertEquals(".GTA", dba.getAlternateAllele(0).getDisplayString());
		dba = new VariantContextDirectedBreakpoint(getContext(), AES(), new VariantContextBuilder(dba).make());
		assertEquals(".GTA", dba.getAlternateAllele(0).getDisplayString());
	}
	@Test
	public void should_set_realignment_failure_flag() {
		VariantContextDirectedBreakpointBuilder builder = new VariantContextDirectedBreakpointBuilder(getContext(), AES())
			.breakend(new BreakendSummary(0, FWD, 1, 1, null), null)
			.realignmentFailed();
		VariantContextDirectedBreakpoint bp = builder.make();
		assertTrue(bp.hasAttribute(VcfAttributes.REALIGNMENT_FAILURE.attribute()));
	}
	@Test
	public void anchor_should_use_reference_base_not_assembly_base() {
		String alt = AB().assemblyBases(B("TTT")).makeVariant().getAlternateAllele(0).getDisplayString();
		assertEquals('A', alt.charAt(0));
	}
	public VariantContextDirectedBreakpoint test_mated_breakend(BreakendDirection direction, boolean realignPositive, String bpString, String realignedCigar, String expectedAllele) {
		VariantContextDirectedBreakpoint dba = AB().direction(direction).anchorLength(0).assemblyBases(B(bpString)).makeVariant();
		dba = new VariantContextDirectedBreakpoint(getContext(), AES(), new VariantContextBuilder(dba).make());
		SAMRecord realigned = Read(1, 10, realignedCigar);
		realigned.setReadBases(realignPositive ? B(bpString) : B(SequenceUtil.reverseComplement(bpString)));
		realigned.setReadNegativeStrandFlag(!realignPositive);
		RealignedBreakpoint rbp = new RealignedBreakpoint(getContext(), dba.getBreakendSummary(), "", realigned);
		VariantContextDirectedBreakpointBuilder builder = new VariantContextDirectedBreakpointBuilder(getContext(), AES(), dba)
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
		VariantContextDirectedBreakpointBuilder builder = new VariantContextDirectedBreakpointBuilder(getContext(), AES())
			.breakend(new BreakendSummary(0, BreakendDirection.Forward, 1, 1, null), null);
		
		assertEquals("BND", builder.make().getAttributeAsString("SVTYPE", ""));
	}
	@Test
	public void breakpoint_should_be_vcf_sv_breakend() {
		VariantContextDirectedBreakpointBuilder builder = new VariantContextDirectedBreakpointBuilder(getContext(), AES())
			.breakpoint(new BreakpointSummary(0, BreakendDirection.Forward, 1, 1, 0, BreakendDirection.Forward, 1, 1, null), null);
		
		assertEquals("BND", builder.make().getAttributeAsString("SVTYPE", ""));
	}
	@Test
	public void evidence_should_round_trip() {
		for (VcfAttributes a : VcfAttributes.evidenceValues()) {
			EvidenceMetrics em = new EvidenceMetrics();
			em.set(a, 1);
			VariantContextDirectedBreakpointBuilder builder = new VariantContextDirectedBreakpointBuilder(getContext(), AES())
				.breakend(new BreakendSummary(0, FWD, 1, 1, null), null)
				.evidence(em);
			builder = new VariantContextDirectedBreakpointBuilder(getContext(), AES(), builder.make());
			
			assertEquals(1, builder.make().getBreakendSummary().evidence.get(a));
		}
	}
	@Test
	public void evidence_score_should_match_vcf_phred_score() {
		EvidenceMetrics m = new EvidenceMetrics(10);
		VariantContextDirectedBreakpointBuilder builder = new VariantContextDirectedBreakpointBuilder(getContext(), AES())
			.breakend(new BreakendSummary(0, FWD, 1, 1, null), null)
			.evidence(m);
		assertEquals(10, builder.make().getPhredScaledQual(), 0);
	}
	@Test
	public void should_round_trip_inexact_breakend() {
		VariantContextDirectedBreakpointBuilder builder = new VariantContextDirectedBreakpointBuilder(getContext(), AES());
		builder.breakend(new BreakendSummary(1, FWD, 2, 4, null), null);
		VariantContextDirectedBreakpoint v = new VariantContextDirectedBreakpointBuilder(getContext(), AES(), builder.make()).make();
		assertEquals(1, v.getBreakendSummary().referenceIndex);
		assertEquals(FWD, v.getBreakendSummary().direction);
		assertEquals(2, v.getBreakendSummary().start);
		assertEquals(4, v.getBreakendSummary().end);
	}
	@Test
	public void should_round_trip_inexact_breakpoint() {
		VariantContextDirectedBreakpointBuilder builder = new VariantContextDirectedBreakpointBuilder(getContext(), AES());
		builder.breakpoint(new BreakpointSummary(1, FWD, 2, 4, 3, BWD, 7, 9, null), null);
		VariantContextDirectedBreakpoint v = new VariantContextDirectedBreakpointBuilder(getContext(), AES(), builder.make()).make();
		BreakpointSummary bp = (BreakpointSummary)v.getBreakendSummary();
		assertEquals(1, bp.referenceIndex);
		assertEquals(FWD, bp.direction);
		assertEquals(2, bp.start);
		assertEquals(4, bp.end);
		assertEquals(3, bp.referenceIndex2);
		assertEquals(BWD, bp.direction2);
		assertEquals(7, bp.start2);
		assertEquals(9, bp.end2);
	}
	@Test
	public void should_write_breakpoint_in_vcf41_mode() {
		ProcessingContext context = new ProcessingContext(
				getFSContext(),
				new ArrayList<Header>(),
				new SoftClipParameters(),
				new AssemblyParameters(),
				new RealignmentParameters(),
				SMALL_FA_FILE, false, true);
		VariantContextDirectedBreakpointBuilder builder = new VariantContextDirectedBreakpointBuilder(context, AES());
		builder.breakend(new BreakendSummary(0,  FWD,  1,  1, null), "ACGT");
		VariantContextDirectedBreakpoint vc = builder.make();
		assertEquals(-1, vc.getAlternateAlleles().get(0).getDisplayString().indexOf("."));
	}
	@Test
	public void vcf41_breakend_should_round_trip() {
		ProcessingContext context = new ProcessingContext(
				getFSContext(),
				new ArrayList<Header>(),
				new SoftClipParameters(),
				new AssemblyParameters(),
				new RealignmentParameters(),
				SMALL_FA_FILE, false, true);
		VariantContextDirectedBreakpointBuilder builder = new VariantContextDirectedBreakpointBuilder(context, AES());
		builder.breakend(new BreakendSummary(0,  FWD,  1,  2, null), "ACGT");
		VariantContextDirectedBreakpoint v = new VariantContextDirectedBreakpointBuilder(context, AES(), builder.make()).make();
		assertEquals(0, v.getBreakendSummary().referenceIndex);
		assertEquals(FWD, v.getBreakendSummary().direction);
		assertEquals(1, v.getBreakendSummary().start);
		assertEquals(2, v.getBreakendSummary().end);
		assertEquals(1, v.getStart());
		assertEquals(1, v.getEnd()); // breakend is called at the first position of the interval
		
		builder = new VariantContextDirectedBreakpointBuilder(context, AES());
		builder.breakend(new BreakendSummary(0,  BWD,  1,  2, null), "ACGT");
		v = new VariantContextDirectedBreakpointBuilder(context, AES(), builder.make()).make();
		assertEquals(0, v.getBreakendSummary().referenceIndex);
		assertEquals(BWD, v.getBreakendSummary().direction);
		assertEquals(1, v.getBreakendSummary().start);
		assertEquals(2, v.getBreakendSummary().end);
		assertEquals(1, v.getStart());
		assertEquals(1, v.getEnd());
	}
	@Test(expected=IllegalStateException.class)
	public void should_require_stop_set() {
		VariantContextDirectedBreakpointBuilder builder = new VariantContextDirectedBreakpointBuilder(getContext(), AES());
		builder.chr("polyA")
			.start(1)
			.alleles("A", "<INS>")
			.make();
	}
}