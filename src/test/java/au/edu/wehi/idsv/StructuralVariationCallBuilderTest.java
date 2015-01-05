package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;
import htsjdk.samtools.SAMRecord;
import htsjdk.variant.vcf.VCFConstants;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import au.edu.wehi.idsv.vcf.VcfAttributes;
import au.edu.wehi.idsv.vcf.VcfFilter;
import au.edu.wehi.idsv.vcf.VcfSvConstants;

import com.google.common.collect.Sets;

public class StructuralVariationCallBuilderTest extends TestHelper {
	public final static BreakpointSummary BP = new BreakpointSummary(0, BWD, 10, 10, 1, BWD, 100, 100);
	public static class sc extends SoftClipEvidence {
		protected sc(int offset, boolean tumour) {
			super(getContext(), SES(tumour), BWD, Read(0, 10, "5S5M"));
			this.offset = offset;
		}
		int offset;
		boolean tumour;
		@Override public int getLocalMapq() { return 1 + offset; }
		@Override public int getLocalBaseLength() { return 2 + offset; }
		@Override public int getLocalMaxBaseQual() { return 4 + offset; }
		@Override public int getLocalTotalBaseQual() { return 5 + offset; }
		@Override public float getBreakendQual() { return 16 + offset; }
		@Override public String getEvidenceID() { return "sc" + Integer.toString(offset); }
	}
	public static class rsc extends RealignedSoftClipEvidence {
		int offset;
		protected rsc(int offset, boolean tumour) {
			super(getContext(), SES(tumour), BWD, Read(0, 10, "5S5M"), onNegative(Read(1, 100, "5M"))[0]);
			this.offset = offset;
		}
		@Override public int getLocalMapq() { return 1 + offset; }
		@Override public int getLocalBaseLength() { return 2 + offset; }
		@Override public int getLocalMaxBaseQual() { return 4 + offset; }
		@Override public int getLocalTotalBaseQual() { return 5 + offset; }
		@Override public int getRemoteMapq() { return 6 + offset; }
		@Override public int getRemoteBaseLength() { return 7 + offset; }
		@Override public int getRemoteBaseCount() { return 8 + offset; }
		@Override public int getRemoteMaxBaseQual() { return 9 + offset; }
		@Override public int getRemoteTotalBaseQual() { return 10 + offset; }
		@Override public float getBreakendQual() { return 111 + offset; }
		@Override public float getBreakpointQual() { return 112 + offset; }
		@Override public String getEvidenceID() { return "rsc" + Integer.toString(offset); }
	}
	public static class rsc_not_supporting_breakpoint extends RealignedSoftClipEvidence {
		int offset;
		protected rsc_not_supporting_breakpoint(int offset, boolean tumour) {
			super(getContext(), SES(tumour), BWD, Read(0, 10, "5S5M"), onNegative(Read(2, 100, "5M"))[0]);
			this.offset = offset;
		}
		@Override public int getLocalMapq() { return 1 + offset; }
		@Override public int getLocalBaseLength() { return 2 + offset; }
		@Override public int getLocalMaxBaseQual() { return 4 + offset; }
		@Override public int getLocalTotalBaseQual() { return 5 + offset; }
		@Override public int getRemoteMapq() { return 6 + offset; }
		@Override public int getRemoteBaseLength() { return 7 + offset; }
		@Override public int getRemoteBaseCount() { return 8 + offset; }
		@Override public int getRemoteMaxBaseQual() { return 9 + offset; }
		@Override public int getRemoteTotalBaseQual() { return 10 + offset; }
		@Override public float getBreakendQual() { return 111 + offset; }
		@Override public float getBreakpointQual() { return 112 + offset; }
		@Override public String getEvidenceID() { return "rsc" + Integer.toString(offset); }
	}
	public static class rrsc extends RealignedRemoteSoftClipEvidence {
		int offset;
		protected rrsc(int offset, boolean tumour) {
			super(getContext(), SES(tumour), BWD, Read(1, 100, "5S5M"), onNegative(Read(0, 10, "5M"))[0]);
			this.offset = offset;
		}
		@Override public int getLocalMapq() { return 1 + offset; }
		@Override public int getLocalBaseLength() { return 2 + offset; }
		@Override public int getLocalMaxBaseQual() { return 4 + offset; }
		@Override public int getLocalTotalBaseQual() { return 5 + offset; }
		@Override public int getRemoteMapq() { return 6 + offset; }
		@Override public int getRemoteBaseLength() { return 7 + offset; }
		@Override public int getRemoteBaseCount() { return 8 + offset; }
		@Override public int getRemoteMaxBaseQual() { return 9 + offset; }
		@Override public int getRemoteTotalBaseQual() { return 10 + offset; }
		@Override public float getBreakendQual() { return 211 + offset; }
		@Override public float getBreakpointQual() { return 212 + offset; }
		@Override public String getEvidenceID() { return "rrsc" + Integer.toString(offset); }
	}
	public static class um extends UnmappedMateReadPair {
		int offset;
		protected um(int offset, boolean tumour) {
			super(OEA(0, 15, "5M", false)[0],
				OEA(0, 15, "5M", false)[1],
				SES(tumour), getContext());
			this.offset = offset;
		}
		@Override public int getLocalMapq() { return 1 + offset; }
		@Override public int getLocalBaseLength() { return 2 + offset; }
		@Override public int getLocalMaxBaseQual() { return 4 + offset; }
		@Override public int getLocalTotalBaseQual() { return 5 + offset; }
		@Override public float getBreakendQual() { return 6 + offset; }
		@Override public String getEvidenceID() { return "um" + Integer.toString(offset); }
	}
	public static class dp extends DiscordantReadPair {
		int offset;
		protected dp(int offset, boolean tumour) {
			super(DP(0, 12, "6M", false, 1, 102, "7M", false)[0],
					DP(0, 12, "6M", false, 1, 102, "7M", false)[1],
					SES(tumour), getContext());
				this.offset = offset;
			}
		@Override public int getLocalMapq() { return 1 + offset; }
		@Override public int getLocalBaseLength() { return 2 + offset; }
		@Override public int getLocalMaxBaseQual() { return 4 + offset; }
		@Override public int getLocalTotalBaseQual() { return 5 + offset; }
		@Override public int getRemoteMapq() { return 6 + offset; }
		@Override public int getRemoteBaseLength() { return 7 + offset; }
		@Override public int getRemoteBaseCount() { return 8 + offset; }
		@Override public int getRemoteMaxBaseQual() { return 9 + offset; }
		@Override public int getRemoteTotalBaseQual() { return 10 + offset; }
		@Override public float getBreakendQual() { return 11 + offset; }
		@Override public float getBreakpointQual() { return 12 + offset; }
		@Override public String getEvidenceID() { return "dp" + Integer.toString(offset); }
	}
	@Test(expected=IllegalArgumentException.class)
	public void should_not_allow_unsupporting_evidence() {
		StructuralVariationCallBuilder cb = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)minimalBreakend()
				.breakend(new BreakendSummary(0, BWD, 1, 10), "").make());
		cb.addEvidence(SCE(BWD, Read(1, 5, "1S2M")));
	}
	@Test
	public void should_set_called_qual() {
		assertTrue(big().hasAttribute(VcfAttributes.CALLED_QUAL.attribute()));
	}
	@Test
	public void should_set_breakend_qual() {
		assertTrue(big().hasAttribute(VcfAttributes.BREAKEND_QUAL.attribute()));
	}
	@Test
	public void should_set_VcfAttribute_REFERENCE_COUNT_READ() {
		assertEquals(7, big().getReferenceReadCount(EvidenceSubset.NORMAL));
		assertEquals(8, big().getReferenceReadCount(EvidenceSubset.TUMOUR));
	}
	@Test
	public void should_set_VcfAttribute_REFERENCE_COUNT_READPAIR() {
		assertEquals(9, big().getReferenceReadPairCount(EvidenceSubset.NORMAL));
		assertEquals(10, big().getReferenceReadPairCount(EvidenceSubset.TUMOUR));
	}
	@Test
	public void should_set_read_pair_evidence_by_breakendpoint_tumournormal() {
		VariantContextDirectedBreakpoint e = (VariantContextDirectedBreakpoint)b(
			new dp(0, true),
			new dp(1, true),
			new dp(2, false),
			new dp(3, false),
			new um(4, true),
			new um(5, true),
			new um(6, true),
			new um(7, false));
		assertEquals(1, e.getBreakendEvidenceCountReadPair(EvidenceSubset.NORMAL));
		assertEquals(3, e.getBreakendEvidenceCountReadPair(EvidenceSubset.TUMOUR));
		assertEquals(2, e.getBreakpointEvidenceCountReadPair(EvidenceSubset.NORMAL));
		assertEquals(2, e.getBreakpointEvidenceCountReadPair(EvidenceSubset.TUMOUR));
	}
	public VariantContextDirectedBreakpoint complex_bp() {
		StructuralVariationCallBuilder cb = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)minimalBreakend()
				.breakpoint(BP, "GT").make());
		cb.addEvidence(new sc(1, true));
		cb.addEvidence(new sc(2, false));
		cb.addEvidence(new sc(3, false));
		cb.addEvidence(new rsc(1, true));
		cb.addEvidence(new rsc(2, true));
		cb.addEvidence(new rsc(3, false));
		cb.addEvidence(new rsc(4, false));
		cb.addEvidence(new rsc(5, false));
		cb.addEvidence(new rrsc(1, true));
		cb.addEvidence(new rrsc(2, true));
		cb.addEvidence(new rrsc(3, true));
		cb.addEvidence(new rrsc(4, false));
		cb.addEvidence(new um(1, false));
		cb.addEvidence(new um(2, false));
		cb.addEvidence(new um(3, false));
		cb.addEvidence(new um(4, false));
		cb.addEvidence(new um(5, true));
		cb.addEvidence(new dp(1, true));
		cb.addEvidence(new dp(2, true));
		cb.addEvidence(AssemblyFactory.incorporateRealignment(getContext(),
				AssemblyFactory.createAnchored(getContext(), AES(), BP.direction, Sets.<DirectedEvidence>newHashSet(
						new rsc(5, false)
						), BP.referenceIndex, BP.end, 1, B("TT"), B("TT"), 1, 2),
				withMapq(44, onNegative(Read(BP.referenceIndex2, BP.start2, "1M")))[0]));
		cb.addEvidence(AssemblyFactory.createAnchored(getContext(), AES(), BP.direction, Sets.<DirectedEvidence>newHashSet(
						new um(6, true)
						), BP.referenceIndex, BP.end, 1, B("TC"), B("TC"), 1, 2));
		cb.addEvidence(AssemblyFactory.incorporateRealignment(getContext(),
				AssemblyFactory.createAnchored(getContext(), AES(), BP.direction, Sets.<DirectedEvidence>newHashSet(
						new rsc(6, false)
						), BP.referenceIndex, BP.end, 1, B("TA"), B("TA"), 1, 2),
				withMapq(1, onNegative(Read(BP.referenceIndex2, BP.start2, "1M")))[0]));
		cb.addEvidence(((RealignedSAMRecordAssemblyEvidence)AssemblyFactory.incorporateRealignment(getContext(),
				AssemblyFactory.createAnchored(getContext(), AES(), BP.direction, Sets.<DirectedEvidence>newHashSet(
						new sc(4, false)), BP.referenceIndex2, BP.end2, 1, B("TT"), B("TT"), 1, 2),
				withMapq(44, onNegative(Read(BP.referenceIndex, BP.start, "1M")))[0])).asRemote());
		return (VariantContextDirectedBreakpoint)cb.make();
	}
	@Test
	public void should_set_BREAKEND_SOFTCLIP() {
		Assert.assertArrayEquals(new int[] { 2,  1, }, (int[])complex_bp().getAttribute(VcfAttributes.BREAKEND_SOFTCLIP_COUNT.attribute()));
		Assert.assertArrayEquals(new float[] { 2*16 + 2+3,  1*16 + 1, }, (float[])complex_bp().getAttribute(VcfAttributes.BREAKEND_SOFTCLIP_QUAL.attribute()), 0);
	}
	@Test
	public void should_set_BREAKPOINT_SOFTCLIP() {
		Assert.assertArrayEquals(new int[] { 3,  2, }, (int[])complex_bp().getAttribute(VcfAttributes.BREAKPOINT_SOFTCLIP_COUNT.attribute()));
		Assert.assertArrayEquals(new float[] { 3*112 + 3+4+5,  2*112 + 1+2, }, (float[])complex_bp().getAttribute(VcfAttributes.BREAKPOINT_SOFTCLIP_QUAL.attribute()), 0);
	}
	@Test
	public void should_set_BREAKPOINT_SOFTCLIP_REMOTE() {
		Assert.assertArrayEquals(new int[] { 1,  3, }, (int[])complex_bp().getAttribute(VcfAttributes.BREAKPOINT_SOFTCLIP_COUNT_REMOTE.attribute()));
		Assert.assertArrayEquals(new float[] { 1*212 + 4, 3*212 + 1+2+3, }, (float[])complex_bp().getAttribute(VcfAttributes.BREAKPOINT_SOFTCLIP_QUAL_REMOTE.attribute()), 0);
	}
	@Test
	public void should_set_BREAKEND_READPAIR() {
		Assert.assertArrayEquals(new int[] { 4,  1, }, (int[])complex_bp().getAttribute(VcfAttributes.BREAKEND_READPAIR_COUNT.attribute()));
		Assert.assertArrayEquals(new float[] { 4*6 + 1+2+3+4,  1*6 + 5, }, (float[])complex_bp().getAttribute(VcfAttributes.BREAKEND_READPAIR_QUAL.attribute()), 0);
	}
	@Test
	public void should_set_BREAKPOINT_READPAIR() {
		Assert.assertArrayEquals(new int[] { 0,  2, }, (int[])complex_bp().getAttribute(VcfAttributes.BREAKPOINT_READPAIR_COUNT.attribute()));
		Assert.assertArrayEquals(new float[] { 0,  2*12 + 1+2, }, (float[])complex_bp().getAttribute(VcfAttributes.BREAKPOINT_READPAIR_QUAL.attribute()), 0);
	}
	@Test
	public void should_not_include_evidence_contributing_to_assembly_in_sc_rp_evidence() {
		should_set_BREAKPOINT_SOFTCLIP_REMOTE();
	}
	@Test
	public void should_set_BREAKEND_ASSEMBLY() {
		Assert.assertEquals(2, (int)(Integer)complex_bp().getAttribute(VcfAttributes.BREAKEND_ASSEMBLY_COUNT.attribute()));
		Assert.assertEquals(14, (float)(Float)complex_bp().getAttribute(VcfAttributes.BREAKEND_ASSEMBLY_QUAL.attribute()), 0);
	}
	@Test
	public void should_set_BREAKPOINT_ASSEMBLY() {
		Assert.assertEquals(1, (int)(Integer)complex_bp().getAttribute(VcfAttributes.BREAKPOINT_ASSEMBLY_COUNT.attribute()));
		Assert.assertEquals(6, (float)(Float)complex_bp().getAttribute(VcfAttributes.BREAKPOINT_ASSEMBLY_QUAL.attribute()), 0);
	}
	@Test
	public void should_set_BREAKPOINT_ASSEMBLY_REMOTE() {
		Assert.assertEquals(1, (int)(Integer)complex_bp().getAttribute(VcfAttributes.BREAKPOINT_ASSEMBLY_COUNT_REMOTE.attribute()));
		Assert.assertEquals(new sc(4, false).getBreakendQual(), (float)(Float)complex_bp().getAttribute(VcfAttributes.BREAKPOINT_ASSEMBLY_QUAL_REMOTE.attribute()), 0);
	}
	/**
	 * We do this to prevent inflation of support for a breakpoint by an assembly that
	 * doesn't not actually support that breakend 
	 */
	@Test
	public void should_exclude_breakpoint_evidence_contributing_to_breakend_assembly_but_not_supporting_breakpoint() {
		StructuralVariationCallBuilder cb = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)minimalBreakend()
				.breakpoint(BP, "GT").make());
		cb.addEvidence(AssemblyFactory.createAnchored(getContext(), AES(), BP.direction, Sets.<DirectedEvidence>newHashSet(
						new rsc_not_supporting_breakpoint(1, true)
						), BP.referenceIndex, BP.end, 1, B("TT"), B("TT"), 1, 2));
		VariantContextDirectedBreakpoint call = (VariantContextDirectedBreakpoint)cb.make();
		Assert.assertArrayEquals(new int[] { 0,  0, }, asIntTN(call.getAttribute(VcfAttributes.BREAKPOINT_SOFTCLIP_COUNT.attribute())));
		Assert.assertArrayEquals(new int[] { 0,  0, }, asIntTN(call.getAttribute(VcfAttributes.BREAKEND_SOFTCLIP_COUNT.attribute())));
	}
	private int[] asIntTN(Object attrValue) {
		if (attrValue == null) return new int[] { 0, 0, };
		return (int[])attrValue;
	}
	@Test
	public void should_count_evidence_once() {
		StructuralVariationCallBuilder cb = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)minimalBreakend()
				.breakpoint(BP, "GT").make());
		cb.addEvidence(new sc(1, true));
		cb.addEvidence(new sc(1, true));
		cb.addEvidence(new rsc(1, true));
		cb.addEvidence(new rsc(1, true));
		cb.addEvidence(new rrsc(1, true));
		cb.addEvidence(new rrsc(1, true));
		cb.addEvidence(new dp(1, true));
		cb.addEvidence(new dp(1, true));
		cb.addEvidence(new um(1, true));
		cb.addEvidence(new um(1, true));
		VariantContextDirectedEvidence v = cb.make();
		Assert.assertArrayEquals(new int[] { 0,  1, }, (int[])v.getAttribute(VcfAttributes.BREAKEND_SOFTCLIP_COUNT.attribute()));
		Assert.assertArrayEquals(new int[] { 0,  1, }, (int[])v.getAttribute(VcfAttributes.BREAKPOINT_SOFTCLIP_COUNT.attribute()));
		Assert.assertArrayEquals(new int[] { 0,  1, }, (int[])v.getAttribute(VcfAttributes.BREAKPOINT_SOFTCLIP_COUNT_REMOTE.attribute()));
		Assert.assertArrayEquals(new int[] { 0,  1, }, (int[])v.getAttribute(VcfAttributes.BREAKEND_READPAIR_COUNT.attribute()));
		Assert.assertArrayEquals(new int[] { 0,  1, }, (int[])v.getAttribute(VcfAttributes.BREAKPOINT_READPAIR_COUNT.attribute()));
	}
	public StructuralVariationCallBuilder cb(DirectedEvidence... evidence) {
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)minimalBreakend()
				.breakpoint(BP, "").make());
		for (DirectedEvidence e : evidence) {
			builder.addEvidence(e);
		}
		return builder;
	}
	public VariantContextDirectedEvidence b(DirectedEvidence... evidence) {
		return cb(evidence).make();
	}
	public VariantContextDirectedEvidence big() {
		ProcessingContext pc = getContext();
		MockSAMEvidenceSource nes = new MockSAMEvidenceSource(pc);
		MockSAMEvidenceSource tes = new MockSAMEvidenceSource(pc);
		tes.isTumour = true;
		StructuralVariationCallBuilder cb = new StructuralVariationCallBuilder(pc, (VariantContextDirectedEvidence)minimalBreakend()
				.breakend(new BreakendSummary(0, BWD, 1, 10), "").make());
		cb.addEvidence(SCE(BWD, nes, Read(0, 5, "1S2M")));
		cb.addEvidence(SCE(BWD, tes, Read(0, 4, "3S4M")));
		cb.addEvidence(SCE(BWD, tes, Read(0, 3, "5S6M")));
		cb.addEvidence(NRRP(nes, DP(0, 11, "5M", false, 1, 100, "6M", true)));
		cb.addEvidence(NRRP(nes, DP(0, 12, "6M", false, 1, 100, "7M", true)));
		cb.addEvidence(NRRP(tes, OEA(0, 13, "8M", false)));
		cb.addEvidence(NRRP(nes, OEA(0, 15, "9M", false)));
		cb.addEvidence(ass1());
		cb.addEvidence(ass2());
		cb.referenceReads(7, 8);
		cb.referenceSpanningPairs(9, 10);
		return cb.make();
	}
	public SAMRecordAssemblyEvidence ass1() {
		ProcessingContext pc = getContext();
		MockSAMEvidenceSource nes = new MockSAMEvidenceSource(pc);
		MockSAMEvidenceSource tes = new MockSAMEvidenceSource(pc);
		tes.isTumour = true;
		Set<DirectedEvidence> support = Sets.newHashSet();
		support.add(SCE(BWD, nes, Read(0, 10, "4S1M")));
		support.add(SCE(BWD, tes, Read(0, 10, "3S5M")));
		support.add(SCE(BWD, tes, Read(0, 10, "3S6M")));
		support.add(NRRP(nes, OEA(0, 15, "5M", false)));
		support.add(NRRP(tes, OEA(0, 16, "5M", false)));
		support.add(NRRP(tes, OEA(0, 17, "5M", false)));
		support.add(NRRP(tes, DP(0, 1, "2M", true, 0, 15, "5M", false)));
		support.add(NRRP(tes, DP(0, 2, "2M", true, 0, 16, "5M", false)));
		support.add(NRRP(nes, DP(0, 3, "2M", true, 0, 17, "10M", false)));
		return AssemblyFactory.createAnchored(pc, AES(), BWD,
				support, 0, 10, 5, B("CGTAAAAT"), new byte[] { 0,1,2,3,4,5,6,7}, 13, 45);
	}
	public RealignedSAMRecordAssemblyEvidence ass2() {
		ProcessingContext pc = getContext();
		MockSAMEvidenceSource nes = new MockSAMEvidenceSource(pc);
		MockSAMEvidenceSource tes = new MockSAMEvidenceSource(pc);
		tes.isTumour = true;
		Set<DirectedEvidence> support = Sets.newHashSet();
		support.add(SCE(BWD, nes, Read(0, 10, "4S1M")));
		support.add(SCE(BWD, tes, Read(0, 10, "3S5M")));
		support.add(SCE(BWD, tes, Read(0, 10, "3S6M")));
		support.add(NRRP(nes, OEA(0, 15, "5M", false)));
		support.add(NRRP(tes, OEA(0, 16, "5M", false)));
		support.add(NRRP(tes, OEA(0, 17, "5M", false)));
		support.add(NRRP(tes, DP(0, 1, "2M", true, 0, 15, "5M", false)));
		support.add(NRRP(tes, DP(0, 2, "2M", true, 0, 16, "5M", false)));
		support.add(NRRP(nes, DP(0, 3, "2M", true, 0, 17, "10M", false)));
		SAMRecordAssemblyEvidence e = AssemblyFactory.createAnchored(pc, AES(), BWD,
				support, 0, 10, 5, B("GGTAAAAC"), new byte[] { 7,6,5,4,3,2,1,0}, 21, 23);
		SAMRecord ra = Read(1, 102, "1S1M1S");
		ra.setReadBases(B("GGT"));
		ra.setMappingQuality(45);
		ra.setBaseQualities(new byte[] { 0,1,2});
		return (RealignedSAMRecordAssemblyEvidence) AssemblyFactory.incorporateRealignment(getContext(), e, ra);
	}
	@Test(expected=IllegalArgumentException.class)
	public void evidence_must_support_call() {
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(),
				(VariantContextDirectedEvidence)minimalBreakend()
					.breakend(new BreakendSummary(0, FWD, 1, 1), null)
					.make());
		builder.addEvidence(SCE(FWD, Read(0, 5, "1M5S")));
		builder.make();
	}
	@Test
	public void should_call_germline_if_evidence_mixed() {
		VariantContextDirectedEvidence e = b(new sc(1, true), new sc(2, false));
		assertFalse(e.hasAttribute(VCFConstants.SOMATIC_KEY));
	}
	@Test
	public void should_call_germline_if_no_normal_evidence() {
		VariantContextDirectedEvidence e = b(new sc(1, true));
		assertFalse(e.hasAttribute(VCFConstants.SOMATIC_KEY));
	}
	@Test
	public void should_somatic_if_tumour_BAF_greater_than_normal() {
		StructuralVariationCallBuilder cb = cb(
				new sc(1, true),
				new sc(2, true),
				new sc(3, true),
				new sc(4, true));
		cb.referenceReads(10, 2);
		cb.referenceSpanningPairs(10, 2);
		assertTrue(cb.make().hasAttribute(VCFConstants.SOMATIC_KEY));
	}
	@Test
	@Ignore("TODO: NYI: Not Yet Implemented")
	public void somatic_p_value_should_be_calculated_from_coverage_at_both_ends_of_the_breakend() {
		fail();
	}
	@Test
	public void should_call_somatic_from_assembly_evidence() {
		ProcessingContext pc = getContext();
		MockSAMEvidenceSource tes = new MockSAMEvidenceSource(pc);
		tes.isTumour = true;
		Set<DirectedEvidence> support = Sets.newHashSet();
		support.add(SCE(BWD, tes, Read(0, 10, "3S5M")));
		support.add(SCE(BWD, tes, Read(0, 10, "3S6M")));
		support.add(NRRP(tes, OEA(0, 16, "5M", false)));
		support.add(NRRP(tes, OEA(0, 17, "5M", false)));
		support.add(NRRP(tes, DP(0, 10, "2M", false, 1, 100, "5M", false)));
		support.add(NRRP(tes, DP(0, 10, "2M", false, 1, 101, "5M", false)));
		AssemblyEvidence ass = AssemblyFactory.createAnchored(pc, AES(), BWD,
				support, 0, 10, 5, B("CGTAAAAT"), new byte[] { 0,1,2,3,4,5,6,7}, 0, 45);
		support.add(ass);
		VariantContextDirectedEvidence e = (VariantContextDirectedEvidence)cb(support.toArray(new DirectedEvidence[0]))
			.referenceReads(10, 10)
			.referenceSpanningPairs(10, 10)
			.make();
		assertTrue(e.hasAttribute(VCFConstants.SOMATIC_KEY));
	}
	@Test
	public void should_filter_small_indels() {
		assertTrue(new StructuralVariationCallBuilder(getContext(),
				(VariantContextDirectedEvidence)minimalBreakend()
					.breakpoint(new BreakpointSummary(new BreakendSummary(0, FWD, 100, 200), new BreakendSummary(0, BWD, 115, 215)), "")
					.make()).make().getFilters().contains(VcfFilter.SMALL_INDEL.filter()));
		assertFalse(new StructuralVariationCallBuilder(getContext(),
				(VariantContextDirectedEvidence)minimalBreakend()
					.breakpoint(new BreakpointSummary(new BreakendSummary(0, FWD, 100, 200), new BreakendSummary(0, BWD, 216, 316)), "")
					.make()).make().getFilters().contains(VcfFilter.SMALL_INDEL.filter()));
	}
	@Test
	public void only_read_pair_should_set_imprecise_header() {
		assertTrue(CallSV(NRRP(OEA(0, 1, "1M", true))).hasAttribute(VcfSvConstants.IMPRECISE_KEY));
		// Should still be inexact event if our CI is 1bp because we don't know if there
		// is any untemplated sequence included in the breakpoint
		assertTrue(CallSV(NRRP(DP(0, 1, "1M", true, 0, 2, "1M", false))).hasAttribute(VcfSvConstants.IMPRECISE_KEY));
	}
	@Test
	public void unanchored_assembly_should_set_imprecise_header() {
		VariantContextDirectedEvidence dba = CallSV(AssemblyFactory.createUnanchored(getContext(), AES(), Sets.<DirectedEvidence>newHashSet(new MockDirectedEvidence(new BreakendSummary(0, FWD, 1, 2))), B("A"), B("A"), 0, 0));
		assertTrue(dba.hasAttribute(VcfSvConstants.IMPRECISE_KEY));
	}
	@Test
	public void should_set_HOMLEN_HOMSEQ_for_microhomology() {
		Assert.fail();
	}
	@Test
	public void should_use_CILEN_tag_for_total_untemplated_sequence() {
		Assert.fail();
	}
	@Test
	public void should_use_INS_UNKNOWN_symbolic_allele_for_unknown_untemplated_sequence() {
		// <INS:UNKNOWN> insertion of unknown sequence and length
		// CIINSLEN confidence interval around insertion  
		Assert.fail();
	}
	@Test
	public void should_use_INS_UNKNOWN_when_bases_missing() {
		// <INS:UNKNOWN> insertion of unknown sequence and length
		// CIINSLEN confidence interval around insertion  
		String localKnownBases = "TT";
		String remoteKnownBases = "CC";
		//
		//                   TODO: work out when to reverse/comp remote untemplated bases
		//                                 vv 
		String localAlt = "ATT<INS:UNKNOWN>CC[remote:100["; 
		String remoteAlt = "ACC<INS:UNKNOWN>TT]local:100]";
		Assert.fail();
	}
	@Test
	public void should_adjust_call_bounds_based_on_best_assembly() {
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
			breakpoint(new BreakpointSummary(0, FWD, 12, 14, 1, BWD, 11, 11), "GTAC");
			phredScore(20);
		}}.make());
		List<DirectedEvidence> evidence = new ArrayList<DirectedEvidence>();
		evidence.add(AssemblyFactory.incorporateRealignment(getContext(),
				AssemblyFactory.createAnchored(getContext(), AES(), FWD, new HashSet<DirectedEvidence>(), 0, 12, 1, B("NTTTT"), B("NTTTT"), 0, 0),
				new SAMRecord(getContext().getBasicSamHeader()) {{
					setReferenceIndex(1);
					setAlignmentStart(11);
					setCigarString("4M");
					setReadNegativeStrandFlag(false);
					setReadBases(B("TTTT"));
					setMappingQuality(21);
					}}));		
		evidence.add(SoftClipEvidence.create(getContext(), SES(), FWD, Read(0, 12, "1M3S"), Read(1, 11, "3M")));
		for (DirectedEvidence e : evidence) {
			builder.addEvidence(e);
		}
		VariantContextDirectedEvidence call = builder.make();
		assertEquals(new BreakpointSummary(0, FWD, 12, 12, 1, BWD, 11, 11), call.getBreakendSummary());
	}
	@Test
	public void should_adjust_call_bounds_based_on_best_soft_clip() {
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
			breakpoint(new BreakpointSummary(0, FWD, 12, 14, 1, BWD, 11, 11), "GTAC");
			phredScore(20);
		}}.make());
		List<DirectedEvidence> evidence = new ArrayList<DirectedEvidence>();
		evidence.add(SoftClipEvidence.create(getContext(), SES(), FWD, Read(0, 12, "1M3S"), withMapq(14, Read(1, 11, "3M"))[0]));
		evidence.add(SoftClipEvidence.create(getContext(), SES(), FWD, Read(0, 13, "1M3S"), withMapq(15, Read(1, 11, "3M"))[0]));
		for (DirectedEvidence e : evidence) {
			builder.addEvidence(e);
		}
		VariantContextDirectedEvidence call = builder.make();
		assertEquals(new BreakpointSummary(0, FWD, 12, 12, 1, BWD, 11, 11), call.getBreakendSummary());
	}
	@Test
	public void should_set_CALLED_QUAL_to_parent_quality_score() {
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
			breakpoint(new BreakpointSummary(0, FWD, 12, 14, 1, BWD, 11, 11), "GTAC");
			phredScore(31);
		}}.make());
		List<DirectedEvidence> evidence = new ArrayList<DirectedEvidence>();
		evidence.add(SoftClipEvidence.create(getContext(), SES(), FWD, Read(0, 12, "1M3S"), withMapq(14, Read(1, 11, "3M"))[0]));
		for (DirectedEvidence e : evidence) {
			builder.addEvidence(e);
		}
		VariantContextDirectedEvidence call = builder.make();
		assertEquals(31, AttributeConverter.asDouble(call.getAttribute(VcfAttributes.CALLED_QUAL.attribute()), 0), 0);
	}
	@Test
	public void should_set_exact_soft_clip_bounds() {
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
			breakpoint(new BreakpointSummary(0, FWD, 11, 12, 0, BWD, 10, 10), "");
			phredScore(10);
		}}.make());
		builder.addEvidence(SoftClipEvidence.create(getContext(), SES(), FWD, withSequence("NNNN", Read(0, 12, "1M3S"))[0], withSequence("NNN", Read(0, 10, "3M"))[0]));
		VariantContextDirectedEvidence de = builder.make();
		assertEquals(new BreakpointSummary(0, FWD, 12, 12, 0, BWD, 10, 10), de.getBreakendSummary());
	}
	@Test(expected=IllegalArgumentException.class)
	public void should_exclude_unsupporting_realigned_soft_clip() {
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
			breakpoint(new BreakpointSummary(0, FWD, 11, 12, 0, BWD, 10, 10), "");
			phredScore(10);
		}}.make());
		builder.addEvidence(SoftClipEvidence.create(getContext(), SES(), FWD, withSequence("NNNN", Read(0, 12, "1M3S"))[0], withSequence("NNN", Read(1, 10, "3M"))[0]));
	}
	@Test(expected=IllegalArgumentException.class)
	public void should_exclude_unsupporting_discordant_read_pair() {
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
			breakpoint(new BreakpointSummary(0, FWD, 11, 12, 0, BWD, 10, 10), "");
			phredScore(10);
		}}.make());
		builder.addEvidence(NRRP(DP(0, 10, "1M", true, 0, 2, "1M", false)));
	}
}