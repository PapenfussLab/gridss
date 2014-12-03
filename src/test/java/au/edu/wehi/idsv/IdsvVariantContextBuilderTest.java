package au.edu.wehi.idsv;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.Header;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import java.util.ArrayList;
import java.util.Collection;

import org.junit.Test;

import au.edu.wehi.idsv.vcf.VcfAttributes;
import au.edu.wehi.idsv.vcf.VcfSvConstants;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

public class IdsvVariantContextBuilderTest extends TestHelper {
	@Test
	public void evidenceID_should_be_ID() {
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(getContext())
			.breakend(new BreakendSummary(0, FWD, 1, 1), null);
		builder.id("testID");
		VariantContextDirectedEvidence bp = (VariantContextDirectedEvidence)builder.make();
		
		assertEquals("testID", bp.getID());
		assertEquals("testID", bp.getEvidenceID());
	}
	@Test
	public void should_lookup_reference_base() {
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(getContext());
		builder.breakend(new BreakendSummary(0, BreakendDirection.Forward, 1, 1), null);
		VariantContextDirectedEvidence bp = (VariantContextDirectedEvidence)builder.make();
		
		assertEquals("A", bp.getReference().getDisplayString());
	}
	@Test
	public void should_round_trip_breakend_quals() {
		byte[] quals = new byte[] { 1, 20, 43 };
		VariantContextDirectedEvidence v = (VariantContextDirectedEvidence)minimalBreakend()
				.breakend(new BreakendSummary(0, FWD, 1, 1), B("AAA"), quals)
				.make();
		v = (VariantContextDirectedEvidence)IdsvVariantContext.create(getContext(), null, v);
		assertArrayEquals(quals, v.getBreakendQuality());
	}
	@Test
	public void should_round_trip_getAssemblyConsensus() {
		VariantContextDirectedEvidence dba = CallSV(AE());
		assertEquals("ATT", dba.getAssemblyConsensus().get(0));
		dba = new VariantContextDirectedEvidence(getContext(), AES(), new VariantContextBuilder(dba).make());
		assertEquals("ATT", dba.getAssemblyConsensus().get(0));
	}
	@Test
	public void getBreakpointSequence_should_get_untemplated_sequence() {
		VariantContextDirectedEvidence dba = CallSV(AE());
		assertEquals("TT", S(dba.getBreakendSequence()));
		dba = new VariantContextDirectedEvidence(getContext(), AES(), new VariantContextBuilder(dba).make());
		assertEquals("TT", S(dba.getBreakendSequence()));
	}
	@Test
	public void should_match_variant_location_f() {
		VariantContextDirectedEvidence dba = CallSV(AssemblyFactory.createAnchored(getContext(), AES(), FWD,
				Sets.<DirectedEvidence>newHashSet(), 0, 10, 1, B("AA"), B("AA"), 0, 0));
		dba = new VariantContextDirectedEvidence(getContext(), AES(), new VariantContextBuilder(dba).make());
		assertEquals("polyA", dba.getChr());
		assertEquals(10, dba.getStart());
		assertEquals(10, dba.getEnd());
	}
	@Test
	public void should_match_variant_location_b() {
		VariantContextDirectedEvidence dba = CallSV(AssemblyFactory.createAnchored(getContext(), AES(), BWD,
				Sets.<DirectedEvidence>newHashSet(), 0, 10, 1, B("AA"), B("AA"), 0, 0));
		dba = new VariantContextDirectedEvidence(getContext(), AES(), new VariantContextBuilder(dba).make());
		assertEquals("polyA", dba.getChr());
		assertEquals(10, dba.getStart());
		assertEquals(10, dba.getEnd());
	}
	@Test
	public void phred_should_be_variant_qual() {
		assertEquals(7.5, minimalBreakend().phredScore(7.5).make().getPhredScaledQual(), 0);
	}
	@Test
	public void should_generate_single_breakend_f() {
		VariantContextDirectedEvidence dba = CallSV(AssemblyFactory.createAnchored(getContext(), AES(), FWD,
				Sets.<DirectedEvidence>newHashSet(), 0, 10, 2, B("NNGT"), B("    "), 0, 0));
		// ref base + breakpoint
		assertEquals("AGT.", dba.getAlternateAllele(0).getDisplayString());
		dba = new VariantContextDirectedEvidence(getContext(), AES(), new VariantContextBuilder(dba).make());
		assertEquals("AGT.", dba.getAlternateAllele(0).getDisplayString());
	}
	@Test
	public void should_generate_single_breakend_b() {
		// ref base + breakpoint
		VariantContextDirectedEvidence dba = CallSV(AssemblyFactory.createAnchored(getContext(), AES(), BWD,
				Sets.<DirectedEvidence>newHashSet(), 0, 10, 2, B("GTNN"), B("    "), 0, 0));
		assertEquals(".GTA", dba.getAlternateAllele(0).getDisplayString());
		dba = new VariantContextDirectedEvidence(getContext(), AES(), new VariantContextBuilder(dba).make());
		assertEquals(".GTA", dba.getAlternateAllele(0).getDisplayString());
	}
	@Test
	public void anchor_should_use_reference_base_not_assembly_base() {
		String alt = CallSV(AssemblyFactory.createAnchored(getContext(), AES(), FWD,
				Sets.<DirectedEvidence>newHashSet(), 0, 10, 1, B("TTT"), B("   "), 0, 0)).getAlternateAllele(0).getDisplayString(); 
		assertEquals('A', alt.charAt(0));
	}
	public VariantContextDirectedEvidence test_mated_breakend(BreakendDirection direction, boolean realignPositive, String bpString, String realignedCigar, String expectedAllele) {
		VariantContextDirectedEvidence dba = CallSV(AssemblyFactory.createUnanchored(getContext(), AES(),
				Sets.<DirectedEvidence>newHashSet(
						NRRP(OEA(0, 1000, "1M", direction == FWD))
						), B(bpString), B(bpString), 0, 0));
				//AB().direction(direction).anchorLength(0).assemblyBases(B(bpString)).makeVariant();
		
		dba = new VariantContextDirectedEvidence(getContext(), AES(), new VariantContextBuilder(dba).make());
		SAMRecord realigned = Read(1, 10, realignedCigar);
		realigned.setReadBases(realignPositive ? B(bpString) : B(SequenceUtil.reverseComplement(bpString)));
		realigned.setReadNegativeStrandFlag(!realignPositive);
		RealignedBreakpoint rbp = new RealignedBreakpoint(getContext(), dba.getBreakendSummary(), "", realigned);
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(getContext(), dba)
			.breakend(rbp.getBreakpointSummary(), rbp.getInsertedSequence());
		dba = (VariantContextDirectedEvidence)builder.make();
		
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
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(getContext())
			.breakend(new BreakendSummary(0, BreakendDirection.Forward, 1, 1), null);
		
		assertEquals("BND", builder.make().getAttributeAsString("SVTYPE", ""));
	}
	@Test
	public void breakpoint_should_be_vcf_sv_breakend() {
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(getContext())
			.breakpoint(new BreakpointSummary(0, BreakendDirection.Forward, 1, 1, 0, BreakendDirection.Forward, 1, 1), null);
		
		assertEquals("BND", builder.make().getAttributeAsString("SVTYPE", ""));
	}
	@Test
	public void should_round_trip_inexact_breakend() {
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(getContext());
		builder.breakend(new BreakendSummary(1, FWD, 2, 4), null);
		VariantContextDirectedEvidence v = (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext(), builder.make()).make();
		assertEquals(1, v.getBreakendSummary().referenceIndex);
		assertEquals(FWD, v.getBreakendSummary().direction);
		assertEquals(2, v.getBreakendSummary().start);
		assertEquals(4, v.getBreakendSummary().end);
	}
	@Test
	public void should_round_trip_inexact_breakpoint() {
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(getContext());
		builder.breakpoint(new BreakpointSummary(1, FWD, 2, 4, 3, BWD, 7, 9), null);
		VariantContextDirectedEvidence v = (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext(), builder.make()).make();
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
	public void should_expose_evidence_source() {
		EvidenceSource source = AES();
		IdsvVariantContextBuilder builder = minimalBreakend();
		builder.source(source);
		assertEquals(source, builder.make().getEvidenceSource());
	}
	@Test
	public void should_write_breakpoint_in_vcf41_mode() {
		ProcessingContext context = new ProcessingContext(
				getFSContext(),
				new ArrayList<Header>(),
				new SoftClipParameters(),
				new ReadPairParameters(),
				new AssemblyParameters(),
				new RealignmentParameters(),
				new VariantCallingParameters(),
				SMALL_FA_FILE, false, true);
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(context);
		builder.breakend(new BreakendSummary(0,  FWD,  1,  1), "ACGT");
		VariantContextDirectedEvidence vc = (VariantContextDirectedEvidence)builder.make();
		assertEquals(-1, vc.getAlternateAlleles().get(0).getDisplayString().indexOf("."));
	}
	@Test
	public void vcf41_breakend_should_round_trip() {
		ProcessingContext context = new ProcessingContext(
				getFSContext(),
				new ArrayList<Header>(),
				new SoftClipParameters(),
				new ReadPairParameters(),
				new AssemblyParameters(),
				new RealignmentParameters(),
				new VariantCallingParameters(),
				SMALL_FA_FILE, false, true);
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(context);
		builder.breakend(new BreakendSummary(0,  FWD,  1,  2), "ACGT");
		VariantContextDirectedEvidence v = (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(context, builder.make()).make();
		assertEquals(0, v.getBreakendSummary().referenceIndex);
		assertEquals(FWD, v.getBreakendSummary().direction);
		assertEquals(1, v.getBreakendSummary().start);
		assertEquals(2, v.getBreakendSummary().end);
		assertEquals(1, v.getStart());
		assertEquals(1, v.getEnd()); // breakend is called at the first position of the interval
		
		builder = new IdsvVariantContextBuilder(context);
		builder.breakend(new BreakendSummary(0,  BWD,  1,  2), "ACGT");
		v = (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(context, builder.make()).make();
		assertEquals(0, v.getBreakendSummary().referenceIndex);
		assertEquals(BWD, v.getBreakendSummary().direction);
		assertEquals(1, v.getBreakendSummary().start);
		assertEquals(2, v.getBreakendSummary().end);
		assertEquals(1, v.getStart());
		assertEquals(1, v.getEnd());
	}
	@Test(expected=IllegalStateException.class)
	public void should_require_stop_set() {
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(getContext());
		builder.chr("polyA")
			.start(1)
			.alleles("A", "<INS>")
			.make();
	}
	@Test
	public void should_not_filter_non_idsv_properties() {
		assertTrue( minimalVariant().attribute("TEST", 0).make().hasAttribute("TEST"));
	}
	@Test
	public void should_not_write_empty_properties() {
		assertFalse( minimalVariant().attribute("LR", null).make().hasAttribute("LR"));
		assertFalse( minimalVariant().attribute("LR", "").make().hasAttribute("LR"));
		assertFalse( minimalVariant().attribute("LR", 0).make().hasAttribute("LR"));
		assertFalse( minimalVariant().attribute("LR", 0L).make().hasAttribute("LR"));
		assertFalse( minimalVariant().attribute("LR", 0f).make().hasAttribute("LR"));
		assertFalse( minimalVariant().attribute("LR", 0d).make().hasAttribute("LR"));
		assertFalse( minimalVariant().attribute("LR", new int[0]).make().hasAttribute("LR"));
		assertFalse( minimalVariant().attribute("LR", new float[0]).make().hasAttribute("LR"));
		assertFalse( minimalVariant().attribute("LR", new double[0]).make().hasAttribute("LR"));
		assertFalse( minimalVariant().attribute("LR", new String[0]).make().hasAttribute("LR"));
		assertFalse( minimalVariant().attribute("LR", Lists.newArrayList()).make().hasAttribute("LR"));
		assertFalse( minimalVariant().attribute("LR", Sets.newHashSet()).make().hasAttribute("LR"));
		
		assertTrue( minimalVariant().attribute("LR", "A").make().hasAttribute("LR"));
		assertTrue( minimalVariant().attribute("LR", new int[] { 1 }).make().hasAttribute("LR"));
		assertTrue( minimalVariant().attribute("LR", new float[] { 1 }).make().hasAttribute("LR"));
		assertTrue( minimalVariant().attribute("LR", new double[] { 1 }).make().hasAttribute("LR"));
		assertTrue( minimalVariant().attribute("LR", new String[] { "s" }).make().hasAttribute("LR"));
		assertTrue( minimalVariant().attribute("LR", Lists.newArrayList("s")).make().hasAttribute("LR"));
		assertTrue( minimalVariant().attribute("LR", Sets.newHashSet("s")).make().hasAttribute("LR"));
	}
	@Test
	public void should_not_write_empty_array_properties() {
		assertFalse( minimalVariant().attribute("LR", new int[] { 0,0 }).make().hasAttribute("LR"));
		assertFalse( minimalVariant().attribute("LR", new float[] { 0,0 }).make().hasAttribute("LR"));
		assertFalse( minimalVariant().attribute("LR", new double[] { 0,0 }).make().hasAttribute("LR"));
		assertFalse( minimalVariant().attribute("LR", new String[] { "", null }).make().hasAttribute("LR"));
		assertFalse( minimalVariant().attribute("LR", new Object[] { 0, 0f, 0d, "", null }).make().hasAttribute("LR"));
		assertFalse( minimalVariant().attribute("LR", Lists.newArrayList("")).make().hasAttribute("LR"));
		assertFalse( minimalVariant().attribute("LR", Lists.<Object>newArrayList(0, 0L, 0f, 0d, "", null)).make().hasAttribute("LR"));
	}
	@SuppressWarnings("unchecked")
	@Test
	public void should_not_trim_array_as_bcf_requires_counts_to_match_header() {
		Object attr = minimalVariant().attribute("LR", ImmutableList.<Integer>of(1,0)).make().getAttribute("LR");
		assertEquals(2, ((Collection<Integer>)attr).size());
	}
	@Test
	public void CIPOS_CIRPOS_should_call_first_position_in_lower_mapped_position_and_corresponding_position_on_remote() {
		// TODO: just call the middle of the interval - that would work very well with SC margins
		VariantContext vc;
		vc = new IdsvVariantContextBuilder(getContext())
			.breakpoint(new BreakpointSummary(0, FWD, 100, 200, 2, BWD, 300, 400), null).make();
		assertEquals(Lists.newArrayList(0,100), AttributeConverter.asIntList(vc.getAttribute(VcfSvConstants.CONFIDENCE_INTERVAL_START_POSITION_KEY)));
		assertEquals(Lists.newArrayList(-100,0), AttributeConverter.asIntList(vc.getAttribute(VcfAttributes.CONFIDENCE_INTERVAL_REMOTE_BREAKEND_START_POSITION_KEY.attribute())));
		
		vc = new IdsvVariantContextBuilder(getContext())
			.breakpoint(new BreakpointSummary(0, FWD, 100, 200, 2, FWD, 300, 400), null).make();
		assertEquals(Lists.newArrayList(0,100), AttributeConverter.asIntList(vc.getAttribute(VcfSvConstants.CONFIDENCE_INTERVAL_START_POSITION_KEY)));
		assertEquals(Lists.newArrayList(0,100), AttributeConverter.asIntList(vc.getAttribute(VcfAttributes.CONFIDENCE_INTERVAL_REMOTE_BREAKEND_START_POSITION_KEY.attribute())));
		
		vc = new IdsvVariantContextBuilder(getContext())
			.breakpoint(new BreakpointSummary(2, BWD, 300, 400, 0, FWD, 100, 200), null).make();
		assertEquals(Lists.newArrayList(-100,0), AttributeConverter.asIntList(vc.getAttribute(VcfSvConstants.CONFIDENCE_INTERVAL_START_POSITION_KEY)));
		assertEquals(Lists.newArrayList(0, 100), AttributeConverter.asIntList(vc.getAttribute(VcfAttributes.CONFIDENCE_INTERVAL_REMOTE_BREAKEND_START_POSITION_KEY.attribute())));
	}
}