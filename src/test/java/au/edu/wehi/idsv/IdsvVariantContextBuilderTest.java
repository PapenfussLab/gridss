package au.edu.wehi.idsv;

import au.edu.wehi.idsv.vcf.VcfInfoAttributes;
import au.edu.wehi.idsv.vcf.VcfSvConstants;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.junit.Ignore;
import org.junit.Test;

import java.util.Collection;

import static org.junit.Assert.*;

public class IdsvVariantContextBuilderTest extends TestHelper {
	@Test
	public void evidenceID_should_be_ID() {
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(getContext())
			.breakend(new BreakendSummary(0, FWD, 1), "");
		builder.id("testID");
		VariantContextDirectedEvidence bp = (VariantContextDirectedEvidence)builder.make();
		
		assertEquals("testID", bp.getID());
		assertEquals("testID", bp.getEvidenceID());
	}
	@Test
	public void should_lookup_reference_base() {
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(getContext());
		builder.breakend(new BreakendSummary(0, BreakendDirection.Forward, 1), "");
		VariantContextDirectedEvidence bp = (VariantContextDirectedEvidence)builder.make();
		
		assertEquals("A", bp.getReference().getDisplayString());
	}
	@Test
	@Ignore("Dropping breakend quals for now")
	public void should_round_trip_breakend_quals() {
		byte[] quals = new byte[] { 1, 20, 43 };
		VariantContextDirectedEvidence v = (VariantContextDirectedEvidence)minimalBreakend()
				.breakend(new BreakendSummary(0, FWD, 1), B("AAA"), quals)
				.make();
		v = (VariantContextDirectedEvidence)IdsvVariantContext.create(getContext().getDictionary(), null, v);
		assertArrayEquals(quals, v.getBreakendQuality());
	}
	@Test
	public void getBreakpointSequence_should_get_untemplated_sequence() {
		VariantContextDirectedEvidence dba = CallSV(AE());
		assertEquals("TT", S(dba.getBreakendSequence()));
		dba = new VariantContextDirectedEvidence(getContext().getDictionary(), AES(), new VariantContextBuilder(dba).make());
		assertEquals("TT", S(dba.getBreakendSequence()));
	}
	@Test
	public void should_match_variant_location_f() {
		VariantContextDirectedEvidence dba = CallSV(AssemblyFactory.createAnchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), FWD,
				null, null, 0, 10, 1, B("AA"), B("AA")));
		dba = new VariantContextDirectedEvidence(getContext().getDictionary(), AES(), new VariantContextBuilder(dba).make());
		assertEquals("polyA", dba.getContig());
		assertEquals(10, dba.getStart());
		assertEquals(10, dba.getEnd());
	}
	@Test
	public void should_match_variant_location_b() {
		VariantContextDirectedEvidence dba = CallSV(AssemblyFactory.createAnchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), BWD,
				null, null, 0, 10, 1, B("AA"), B("AA")));
		dba = new VariantContextDirectedEvidence(getContext().getDictionary(), AES(), new VariantContextBuilder(dba).make());
		assertEquals("polyA", dba.getContig());
		assertEquals(10, dba.getStart());
		assertEquals(10, dba.getEnd());
	}
	@Test
	public void phred_should_be_variant_qual() {
		assertEquals(7.5, minimalBreakend().phredScore(7.5).make().getPhredScaledQual(), 0);
	}
	@Test
	public void should_generate_single_breakend_f() {
		VariantContextDirectedEvidence dba = CallSV(AssemblyFactory.createAnchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), FWD,
				null, null, 0, 10, 2, B("NNGT"), B("    ")));
		// ref base + breakpoint
		assertEquals("AGT.", dba.getAlternateAllele(0).getDisplayString());
		dba = new VariantContextDirectedEvidence(getContext().getDictionary(), AES(), new VariantContextBuilder(dba).make());
		assertEquals("AGT.", dba.getAlternateAllele(0).getDisplayString());
	}
	@Test
	public void should_generate_single_breakend_b() {
		// ref base + breakpoint
		VariantContextDirectedEvidence dba = CallSV(AssemblyFactory.createAnchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), BWD,
				null, null, 0, 10, 2, B("GTNN"), B("    ")));
		assertEquals(".GTA", dba.getAlternateAllele(0).getDisplayString());
		dba = new VariantContextDirectedEvidence(getContext().getDictionary(), AES(), new VariantContextBuilder(dba).make());
		assertEquals(".GTA", dba.getAlternateAllele(0).getDisplayString());
	}
	@Test
	public void anchor_should_use_reference_base_not_assembly_base() {
		String alt = CallSV(AssemblyFactory.createAnchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), FWD,
				null, null, 0, 10, 1, B("TTT"), B("   "))).getAlternateAllele(0).getDisplayString();
		assertEquals('A', alt.charAt(0));
	}
	/*// These are actually RealignedBreakpoint tests
	public VariantContextDirectedEvidence test_mated_breakend(BreakendDirection direction, boolean realignPositive, String bpString, String realignedCigar, String expectedAllele) {
		SAMRecord realigned = Read(1, 10, realignedCigar);
		realigned.setReadBases(realignPositive ? B(bpString) : B(SequenceUtil.reverseComplement(bpString)));
		realigned.setReadNegativeStrandFlag(!realignPositive);
		int fragSize = realigned.getCigar().getReferenceLength() - 1; // temp hack until better upstream calculation is performed
		VariantContextDirectedEvidence dba = CallSV(AssemblyFactory.createUnanchoredBreakend(getContext(), AES(fragSize),
				new BreakendSummary(0, direction, 1000), null, //NRRP(SES(fragSize), OEA(0, 1000, "1M", direction == FWD))
				B(bpString), B(bpString), new int[] {0, 0}));
				//AB().direction(direction).anchorLength(0).assemblyBases(B(bpString)).makeVariant();
		
		dba = new VariantContextDirectedEvidence(getContext(), AES(fragSize), new VariantContextBuilder(dba).make());
		
		RealignedBreakpoint rbp = RealignedBreakpoint.create(getContext().getReference(), getContext().getDictionary(), dba.getBreakendSummary(), "", realigned);
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
	*/
	@Test
	public void breakend_should_be_vcf_sv_breakend() {
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(getContext())
			.breakend(new BreakendSummary(0, BreakendDirection.Forward, 1), "");
		
		assertEquals("BND", builder.make().getAttributeAsString("SVTYPE", ""));
	}
	@Test
	public void breakpoint_should_be_vcf_sv_breakend() {
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(getContext())
			.breakpoint(new BreakpointSummary(0, BreakendDirection.Forward, 1, 0, BreakendDirection.Forward, 1), "");
		
		assertEquals("BND", builder.make().getAttributeAsString("SVTYPE", ""));
	}
	@Test
	public void should_round_trip_inexact_breakend() {
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(getContext());
		builder.breakend(new BreakendSummary(1, FWD, 3, 2, 7), "");
		VariantContextDirectedEvidence v = (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext(), builder.make()).make();
		assertEquals(1, v.getBreakendSummary().referenceIndex);
		assertEquals(FWD, v.getBreakendSummary().direction);
		assertEquals(2, v.getBreakendSummary().start);
		assertEquals(7, v.getBreakendSummary().end);
		assertEquals(3, v.getBreakendSummary().nominal);
	}
	@Test
	public void should_round_trip_inexact_breakpoint() {
		for (final BreakpointSummary bp : new BreakpointSummary[] {
			new BreakpointSummary(0, FWD, 12, 11, 12, 0, BWD, 8, 5, 15),
			new BreakpointSummary(0, FWD, 10, 10, 10, 0, BWD, 11, 10, 15),
			new BreakpointSummary(0, FWD, 11, 10, 12, 0, BWD, 10, 10, 10),
			new BreakpointSummary(0, FWD, 12, 11, 12, 0, FWD, 13, 10, 20),
		}) {
			assertEquals(bp, ((VariantContextDirectedEvidence)(new IdsvVariantContextBuilder(getContext()) {{
				breakpoint(bp, "GTAC");
			}}.make())).getBreakendSummary());
		}
	}
	@Test
	public void should_expose_evidence_source() {
		EvidenceSource source = AES();
		IdsvVariantContextBuilder builder = minimalBreakend();
		builder.source(source);
		assertEquals(source, builder.make().getEvidenceSource());
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
		assertTrue( minimalVariant().attribute("UNUSED_FIELD", 0).make().hasAttribute("UNUSED_FIELD"));
	}
	@Test
	public void should_write_empty_properties() {
		assertEquals(gridss.Defaults.WRITE_ZERO_OR_EMTPY_VCF_FIELDS, minimalVariant().attribute("TEST", null).make().hasAttribute("TEST"));
		assertEquals(gridss.Defaults.WRITE_ZERO_OR_EMTPY_VCF_FIELDS, minimalVariant().attribute("TEST", "").make().hasAttribute("TEST"));
		assertEquals(gridss.Defaults.WRITE_ZERO_OR_EMTPY_VCF_FIELDS, minimalVariant().attribute("TEST", 0).make().hasAttribute("TEST"));
		assertEquals(gridss.Defaults.WRITE_ZERO_OR_EMTPY_VCF_FIELDS, minimalVariant().attribute("TEST", 0L).make().hasAttribute("TEST"));
		assertEquals(gridss.Defaults.WRITE_ZERO_OR_EMTPY_VCF_FIELDS, minimalVariant().attribute("TEST", 0f).make().hasAttribute("TEST"));
		assertEquals(gridss.Defaults.WRITE_ZERO_OR_EMTPY_VCF_FIELDS, minimalVariant().attribute("TEST", 0d).make().hasAttribute("TEST"));
		assertEquals(gridss.Defaults.WRITE_ZERO_OR_EMTPY_VCF_FIELDS, minimalVariant().attribute("TEST", new int[0]).make().hasAttribute("TEST"));
		assertEquals(gridss.Defaults.WRITE_ZERO_OR_EMTPY_VCF_FIELDS, minimalVariant().attribute("TEST", new float[0]).make().hasAttribute("TEST"));
		assertEquals(gridss.Defaults.WRITE_ZERO_OR_EMTPY_VCF_FIELDS, minimalVariant().attribute("TEST", new double[0]).make().hasAttribute("TEST"));
		assertEquals(gridss.Defaults.WRITE_ZERO_OR_EMTPY_VCF_FIELDS, minimalVariant().attribute("TEST", new String[0]).make().hasAttribute("TEST"));
		assertEquals(gridss.Defaults.WRITE_ZERO_OR_EMTPY_VCF_FIELDS, minimalVariant().attribute("TEST", Lists.newArrayList()).make().hasAttribute("TEST"));
		assertEquals(gridss.Defaults.WRITE_ZERO_OR_EMTPY_VCF_FIELDS, minimalVariant().attribute("TEST", Sets.newHashSet()).make().hasAttribute("TEST"));
		
		assertTrue( minimalVariant().attribute("TEST", "A").make().hasAttribute("TEST"));
		assertTrue( minimalVariant().attribute("TEST", new int[] { 1 }).make().hasAttribute("TEST"));
		assertTrue( minimalVariant().attribute("TEST", new float[] { 1 }).make().hasAttribute("TEST"));
		assertTrue( minimalVariant().attribute("TEST", new double[] { 1 }).make().hasAttribute("TEST"));
		assertTrue( minimalVariant().attribute("TEST", new String[] { "s" }).make().hasAttribute("TEST"));
		assertTrue( minimalVariant().attribute("TEST", Lists.newArrayList("s")).make().hasAttribute("TEST"));
		assertTrue( minimalVariant().attribute("TEST", Sets.newHashSet("s")).make().hasAttribute("TEST"));
	}
	@Test
	public void should_write_empty_array_properties() {
		assertEquals(gridss.Defaults.WRITE_ZERO_OR_EMTPY_VCF_FIELDS, minimalVariant().attribute("TEST", new int[] { 0,0 }).make().hasAttribute("TEST"));
		assertEquals(gridss.Defaults.WRITE_ZERO_OR_EMTPY_VCF_FIELDS, minimalVariant().attribute("TEST", new float[] { 0,0 }).make().hasAttribute("TEST"));
		assertEquals(gridss.Defaults.WRITE_ZERO_OR_EMTPY_VCF_FIELDS, minimalVariant().attribute("TEST", new double[] { 0,0 }).make().hasAttribute("TEST"));
		assertEquals(gridss.Defaults.WRITE_ZERO_OR_EMTPY_VCF_FIELDS, minimalVariant().attribute("TEST", new String[] { "", null }).make().hasAttribute("TEST"));
		assertEquals(gridss.Defaults.WRITE_ZERO_OR_EMTPY_VCF_FIELDS, minimalVariant().attribute("TEST", new Object[] { 0, 0f, 0d, "", null }).make().hasAttribute("TEST"));
		assertEquals(gridss.Defaults.WRITE_ZERO_OR_EMTPY_VCF_FIELDS, minimalVariant().attribute("TEST", Lists.newArrayList("")).make().hasAttribute("TEST"));
		assertEquals(gridss.Defaults.WRITE_ZERO_OR_EMTPY_VCF_FIELDS, minimalVariant().attribute("TEST", Lists.<Object>newArrayList(0, 0L, 0f, 0d, "", null)).make().hasAttribute("TEST"));
	}
	@SuppressWarnings("unchecked")
	@Test
	public void should_not_trim_array_as_bcf_requires_counts_to_match_header() {
		Object attr = minimalVariant().attribute("LR", ImmutableList.<Integer>of(1,0)).make().getAttribute("LR");
		assertEquals(2, ((Collection<Integer>)attr).size());
	}
	@Test
	public void CIPOS_CIRPOS_should_call_nominal_position() {
		VariantContext vc;
		vc = new IdsvVariantContextBuilder(getContext())
			.breakpoint(new BreakpointSummary(0, FWD, 101, 100, 102, 2, BWD, 301, 300, 302), "").make();
		assertEquals(Lists.newArrayList(-1,1), AttributeConverter.asIntList(vc.getAttribute(VcfSvConstants.CONFIDENCE_INTERVAL_START_POSITION_KEY)));
		assertEquals(Lists.newArrayList(-1,1), AttributeConverter.asIntList(vc.getAttribute(VcfInfoAttributes.CONFIDENCE_INTERVAL_REMOTE_BREAKEND_START_POSITION_KEY.attribute())));
		
		vc = new IdsvVariantContextBuilder(getContext())
			.breakpoint(new BreakpointSummary(0, FWD, 150, 100, 200, 2, FWD, 350, 300, 400), "").make();
		assertEquals(Lists.newArrayList(-50,50), AttributeConverter.asIntList(vc.getAttribute(VcfSvConstants.CONFIDENCE_INTERVAL_START_POSITION_KEY)));
		assertEquals(Lists.newArrayList(-50,50), AttributeConverter.asIntList(vc.getAttribute(VcfInfoAttributes.CONFIDENCE_INTERVAL_REMOTE_BREAKEND_START_POSITION_KEY.attribute())));
		
		vc = new IdsvVariantContextBuilder(getContext())
			.breakpoint(new BreakpointSummary(2, BWD, 350, 300, 400, 0, FWD, 150, 100, 200), "").make();
		assertEquals(Lists.newArrayList(-50,50), AttributeConverter.asIntList(vc.getAttribute(VcfSvConstants.CONFIDENCE_INTERVAL_START_POSITION_KEY)));
		assertEquals(Lists.newArrayList(-50,50), AttributeConverter.asIntList(vc.getAttribute(VcfInfoAttributes.CONFIDENCE_INTERVAL_REMOTE_BREAKEND_START_POSITION_KEY.attribute())));
		
		vc = new IdsvVariantContextBuilder(getContext())
				.breakpoint(new BreakpointSummary(0, FWD, 101, 100, 105, 2, BWD, 300, 300, 305), "").make();
			assertEquals(Lists.newArrayList(-1,4), AttributeConverter.asIntList(vc.getAttribute(VcfSvConstants.CONFIDENCE_INTERVAL_START_POSITION_KEY)));
			assertEquals(Lists.newArrayList(0,5), AttributeConverter.asIntList(vc.getAttribute(VcfInfoAttributes.CONFIDENCE_INTERVAL_REMOTE_BREAKEND_START_POSITION_KEY.attribute())));
	}
}