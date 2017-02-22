package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;

import au.edu.wehi.idsv.sam.ChimericAlignment;
import au.edu.wehi.idsv.vcf.VcfFormatAttributes;
import au.edu.wehi.idsv.vcf.VcfInfoAttributes;
import au.edu.wehi.idsv.vcf.VcfSvConstants;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.variant.vcf.VCFConstants;

public class StructuralVariationCallBuilderTest extends TestHelper {
	public final static BreakpointSummary BP = new BreakpointSummary(0, BWD, 10, 1, BWD, 100);
	public static class sc extends SoftClipEvidence {
		protected sc(int offset, boolean tumour) {
			super(SES(tumour), withSequence("NNNNNNNNNN", Read(0, 10, "5S5M"))[0],
				new BreakendSummary(0, BWD, 10, 10, 10),
				0, 5, 5, 10, 0);
			this.offset = offset;
		}
		int offset;
		boolean tumour;
		@Override public int getLocalMapq() { return 1 + offset; }
		@Override public float getBreakendQual() { return 16 + offset; }
		@Override public String getEvidenceID() { return "sc" + Integer.toString(offset); }
	}
	public static class rsc extends SplitReadEvidence {
		int offset;
		protected rsc(int offset, boolean tumour) {
			super(SES(tumour), withAttr("SA", "polyACGT,100,5M,-,0,0", withSequence("NNNNNNNNNN", Read(0, 10, "5S5M")))[0],
					new BreakpointSummary(0, BWD, 10, 1, BWD, 100),
					0, 5, 5, 5, 5, 10, new ChimericAlignment("polyACGT,100,-,5M,0,0"),0, 0);
			this.offset = offset;
		}
		@Override public int getLocalMapq() { return 1 + offset; }
		@Override public int getRemoteMapq() { return 6 + offset; }
		@Override public float getBreakendQual() { return 111 + offset; }
		@Override public float getBreakpointQual() { return 112 + offset; }
		@Override public String getEvidenceID() { return "rsc" + Integer.toString(offset); }
	}
	public static SAMRecord[] asSupplementary(SAMRecord... r) {
		for (SAMRecord rec : r) {
			rec.setSupplementaryAlignmentFlag(true);
		}
		return r;
	}
	public static class rsc_not_supporting_breakpoint extends SplitReadEvidence {
		int offset;
		protected rsc_not_supporting_breakpoint(int offset, boolean tumour) {
			super(SES(tumour), withAttr("SA", "random,100,5M,-,0,0", withSequence("NNNNNNNNNN", Read(0, 10, "5S5M")))[0],
					new BreakpointSummary(0, BWD, 10, 2, BWD, 100),
					0, 5, 5, 5, 5, 10, new ChimericAlignment("random,100,-,5M,0,0"),0, 0);
		}
		@Override public int getLocalMapq() { return 1 + offset; }
		@Override public int getRemoteMapq() { return 6 + offset; }
		@Override public float getBreakendQual() { return 111 + offset; }
		@Override public float getBreakpointQual() { return 112 + offset; }
		@Override public String getEvidenceID() { return "rsc" + Integer.toString(offset); }
	}
	public static class rrsc extends SplitReadEvidence {
		int offset;
		protected rrsc(int offset, boolean tumour) {
			super(SES(tumour), asSupplementary(withAttr("SA", "polyACGT,100,-,5M,0,0", withSequence("NNNNNNNNNN", Read(0, 10, "5S5M"))))[0],
					new BreakpointSummary(0, BWD, 10, 1, BWD, 100),
					0, 5, 5, 5, 5, 10, new ChimericAlignment("polyACGT,100,-,5M,0,0"),0, 0);
			this.offset = offset;
		}
		@Override public int getLocalMapq() { return 1 + offset; }
		@Override public int getRemoteMapq() { return 6 + offset; }
		@Override public float getBreakendQual() { return 211 + offset; }
		@Override public float getBreakpointQual() { return 212 + offset; }
		@Override public String getEvidenceID() { return "Rrsc" + Integer.toString(offset); }
	}
	public static class um extends UnmappedMateReadPair {
		int offset;
		protected um(int offset, boolean tumour) {
			super(OEA(0, 15, "5M", false)[0],
				OEA(0, 15, "5M", false)[1],
				SES(tumour));
			this.offset = offset;
		}
		@Override public int getLocalMapq() { return 1 + offset; }
		@Override public float getBreakendQual() { return 6 + offset; }
		@Override public String getEvidenceID() { return "um" + Integer.toString(offset); }
	}
	public static class dp extends DiscordantReadPair {
		int offset;
		protected dp(int offset, boolean tumour) {
			super(DP(0, 12, "6M", false, 1, 102, "7M", false)[0],
					DP(0, 12, "6M", false, 1, 102, "7M", false)[1],
					SES(tumour));
				this.offset = offset;
			}
		@Override public int getLocalMapq() { return 1 + offset; }
		@Override public int getRemoteMapq() { return 6 + offset; }
		@Override public float getBreakendQual() { return 11 + offset; }
		@Override public float getBreakpointQual() { return 12 + offset; }
		@Override public String getEvidenceID() { return "dp" + Integer.toString(offset); }
	}
	@Test(expected=IllegalArgumentException.class)
	public void should_not_allow_unsupporting_evidence() {
		StructuralVariationCallBuilder cb = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)minimalBreakend()
				.breakend(new BreakendSummary(0, BWD, 1, 1, 10), "").make());
		cb.addEvidence(SCE(BWD, Read(1, 5, "1S2M")));
	}
	@Test
	public void should_set_called_qual() {
		assertTrue(big().hasAttribute(VcfInfoAttributes.CALLED_QUAL.attribute()));
	}
	@Test
	public void should_set_breakend_qual() {
		assertTrue(big().hasAttribute(VcfInfoAttributes.BREAKEND_QUAL.attribute()));
	}
	@Test
	public void should_set_VcfAttribute_REFERENCE_COUNT_READ() {
		assertEquals(7, big().getReferenceReadCount(0));
		assertEquals(8, big().getReferenceReadCount(1));
	}
	@Test
	public void should_set_VcfAttribute_REFERENCE_COUNT_READPAIR() {
		assertEquals(9, big().getReferenceReadPairCount(0));
		assertEquals(10, big().getReferenceReadPairCount(1));
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
		assertEquals(1, e.getBreakendEvidenceCountReadPair(0));
		assertEquals(3, e.getBreakendEvidenceCountReadPair(1));
		assertEquals(2, e.getBreakpointEvidenceCountReadPair(0));
		assertEquals(2, e.getBreakpointEvidenceCountReadPair(1));
	}
	public VariantContextDirectedBreakpoint complex_bp() {
		ProcessingContext pc = getContext();
		AssemblyEvidenceSource aes = AES(pc);
		StructuralVariationCallBuilder cb = new StructuralVariationCallBuilder(pc, (VariantContextDirectedEvidence)minimalBreakend()
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
		
		List<DirectedEvidence> support = Lists.<DirectedEvidence>newArrayList(new rsc(5, false)); 
		SAMRecord ass = AssemblyFactory.createAnchoredBreakend(pc, aes, new SequentialIdGenerator("asm"), BP.direction, support, BP.referenceIndex, BP.end, 1, B("TT"), B("TT"));
		cb.addEvidence(incorporateRealignment(AES(), ass, ImmutableList.of(withMapq(44, onNegative(Read(BP.referenceIndex2, BP.start2, "1M")))[0])));
		
		support = Lists.<DirectedEvidence>newArrayList(new um(6, true)); 
		ass = AssemblyFactory.createAnchoredBreakend(pc, aes, new SequentialIdGenerator("asm2"), BP.direction, support, BP.referenceIndex, BP.end, 1, B("TT"), B("TT"));
		cb.addEvidence(asAssemblyEvidence(ass));
		
		support = Lists.<DirectedEvidence>newArrayList(new rsc(6, false));
		ass = AssemblyFactory.createAnchoredBreakend(pc, aes, new SequentialIdGenerator("asm3"), BP.direction, support, BP.referenceIndex, BP.end, 1, B("TA"), B("TA"));
		cb.addEvidence(incorporateRealignment(aes, ass, ImmutableList.of(withMapq(1, onNegative(Read(BP.referenceIndex2, BP.start2, "1M")))[0])));
		
		support = Lists.<DirectedEvidence>newArrayList(new sc(4, false));
		ass = AssemblyFactory.createAnchoredBreakend(pc, aes, new SequentialIdGenerator("asm4"), BP.direction, support, BP.referenceIndex2, BP.end2, 1, B("TT"), B("TT"));
		SAMRecord remote = withMapq(44, onNegative(Read(BP.referenceIndex, BP.start, "1M")))[0];
		incorporateRealignment(aes, ass, ImmutableList.of(remote));
		cb.addEvidence(asAssemblyEvidence(remote));
		
		return (VariantContextDirectedBreakpoint)cb.make();
	}
	private void assertAttr(VcfInfoAttributes infoAttr, VcfFormatAttributes formatAttr, int[] expected, VariantContextDirectedEvidence bp) {
		Assert.assertEquals(Arrays.stream(expected).sum(), bp.getAttribute(infoAttr.attribute()));
		for (int i = 0; i < expected.length; i++) {
			Assert.assertEquals((Integer)expected[i], bp.getGenotype(i).getExtendedAttribute(formatAttr.attribute(), Integer.MIN_VALUE + 1));
		}
	}
	private void assertAttr(VcfInfoAttributes infoAttr, VcfFormatAttributes formatAttr, double[] expected, VariantContextDirectedEvidence bp) {
		Assert.assertEquals(Arrays.stream(expected).sum(), bp.getAttribute(infoAttr.attribute()));
		for (int i = 0; i < expected.length; i++) {
			Assert.assertEquals((Double)expected[i], bp.getGenotype(i).getExtendedAttribute(formatAttr.attribute(), Integer.MIN_VALUE + 1));
		}
	}
	@Test
	public void should_set_BREAKEND_SOFTCLIP() {
		assertAttr(VcfInfoAttributes.BREAKEND_SOFTCLIP_COUNT, VcfFormatAttributes.BREAKEND_SOFTCLIP_COUNT, new int[] { 2,  1, }, complex_bp());
		assertAttr(VcfInfoAttributes.BREAKEND_SOFTCLIP_QUAL, VcfFormatAttributes.BREAKEND_SOFTCLIP_QUAL, new double[] { 2*16 + 2+3,  1*16 + 1, }, complex_bp());
	}
	@Test
	public void should_set_BREAKPOINT_SPLITREAD() {
		assertAttr(VcfInfoAttributes.BREAKPOINT_SPLITREAD_COUNT, VcfFormatAttributes.BREAKPOINT_SPLITREAD_COUNT, new int[] { 4,  5, }, complex_bp());
		assertAttr(VcfInfoAttributes.BREAKPOINT_SPLITREAD_QUAL, VcfFormatAttributes.BREAKPOINT_SPLITREAD_QUAL, new double[] { 3*112 + 3+4+5 + 212+4,  2*112 + 1+2 + 3*212+1+2+3, }, complex_bp());
	}
	@Test
	public void should_set_BREAKEND_READPAIR() {
		assertAttr(VcfInfoAttributes.BREAKEND_UNMAPPEDMATE_COUNT, VcfFormatAttributes.BREAKEND_UNMAPPEDMATE_COUNT, new int[] { 4,  1, }, complex_bp());
		assertAttr(VcfInfoAttributes.BREAKEND_UNMAPPEDMATE_QUAL, VcfFormatAttributes.BREAKEND_UNMAPPEDMATE_QUAL, new double[]  { 4*6 + 1+2+3+4,  1*6 + 5, }, complex_bp());
	}
	@Test
	public void should_set_BREAKPOINT_READPAIR() {
		assertAttr(VcfInfoAttributes.BREAKPOINT_READPAIR_COUNT, VcfFormatAttributes.BREAKPOINT_READPAIR_COUNT, new int[] { 0,  2, }, complex_bp());
		assertAttr(VcfInfoAttributes.BREAKPOINT_READPAIR_QUAL, VcfFormatAttributes.BREAKPOINT_READPAIR_QUAL, new double[]  { 0,  2*12 + 1+2, }, complex_bp());
	}
	@Test
	public void should_set_BREAKEND_ASSEMBLY() {
		Assert.assertEquals(1, (int)(Integer)complex_bp().getAttribute(VcfInfoAttributes.BREAKEND_ASSEMBLY_COUNT.attribute()));
		Assert.assertEquals(7, (double)(Double)complex_bp().getAttribute(VcfInfoAttributes.BREAKEND_ASSEMBLY_QUAL.attribute()), 0);
	}
	@Test
	public void should_set_BREAKPOINT_ASSEMBLY() {
		Assert.assertEquals(2, (int)(Integer)complex_bp().getAttribute(VcfInfoAttributes.BREAKPOINT_ASSEMBLY_COUNT.attribute()));
		Assert.assertEquals(7, (double)(Double)complex_bp().getAttribute(VcfInfoAttributes.BREAKPOINT_ASSEMBLY_QUAL.attribute()), 0);
	}
	@Test
	public void should_set_BREAKPOINT_ASSEMBLY_REMOTE() {
		Assert.assertEquals(1, (int)(Integer)complex_bp().getAttribute(VcfInfoAttributes.BREAKPOINT_ASSEMBLY_COUNT_REMOTE.attribute()));
		Assert.assertEquals(5, (double)(Double)complex_bp().getAttribute(VcfInfoAttributes.BREAKPOINT_ASSEMBLY_QUAL_REMOTE.attribute()), 0);
	}
	@Test
	public void should_set_BREAKPOINT_INDEL() {
		IndelEvidence ie = IE(SES(true), Read(0, 1, "10M10D10M"));
		StructuralVariationCallBuilder cb = new StructuralVariationCallBuilder(getContext(), BP("variant", new BreakpointSummary(0, FWD, 10, 0, BWD, 21)));
		cb.addEvidence(ie);
		VariantContextDirectedEvidence vc = cb.make();
		
		assertAttr(VcfInfoAttributes.BREAKPOINT_INDEL_COUNT, VcfFormatAttributes.BREAKPOINT_INDEL_COUNT, new int[] { 0,  1, }, vc);
		assertAttr(VcfInfoAttributes.BREAKPOINT_INDEL_QUAL, VcfFormatAttributes.BREAKPOINT_INDEL_QUAL, new double[]  { 0,  ie.getBreakpointQual(), }, vc);
	}
	@Test
	public void should_set_VcfAttribute_BREAKPOINT_ASSEMBLY_READPAIR_COUNT() {
		ProcessingContext pc = getContext();
		AssemblyEvidenceSource aes = AES(pc);
		StructuralVariationCallBuilder cb = new StructuralVariationCallBuilder(pc, (VariantContextDirectedEvidence)minimalBreakend()
				.breakpoint(BP, "GT").make());
		cb.addEvidence(new dp(1, true));
		cb.addEvidence(new dp(2, true));
		cb.addEvidence(new dp(3, false));
		
		List<DirectedEvidence> support = Lists.<DirectedEvidence>newArrayList(
				new dp(2, true),
				new dp(4, true).asRemote()); 
		SAMRecord ass = AssemblyFactory.createAnchoredBreakend(pc, aes, new SequentialIdGenerator("asm"), BP.direction, support, BP.referenceIndex, BP.end, 1, B("TT"), B("TT"));
		cb.addEvidence(incorporateRealignment(AES(), ass, ImmutableList.of(withMapq(44, onNegative(Read(BP.referenceIndex2, BP.start2, "1M")))[0])));
	
		support = Lists.<DirectedEvidence>newArrayList(
				// none support since this is the remote breakend
				new dp(3, false),
				new dp(5, true),
				new dp(6, true));
		ass = AssemblyFactory.createAnchoredBreakend(pc, aes, new SequentialIdGenerator("asm"), BP.direction, support, BP.referenceIndex2, BP.end2, 1, B("TT"), B("TT"));
		SAMRecord remote = withMapq(44, onNegative(Read(BP.referenceIndex, BP.start, "1M")))[0];
		incorporateRealignment(aes, ass, ImmutableList.of(remote));
		cb.addEvidence(asAssemblyEvidence(remote));
		
		VariantContextDirectedBreakpoint bp = (VariantContextDirectedBreakpoint)cb.make();
		
		assertAttr(VcfInfoAttributes.BREAKPOINT_ASSEMBLY_READPAIR_COUNT, VcfFormatAttributes.BREAKPOINT_ASSEMBLY_READPAIR_COUNT, new int[] { 1, 4, }, bp);
	}
	@Test
	public void should_set_VcfAttribute_BREAKPOINT_ASSEMBLY_READ_COUNT() {
		ProcessingContext pc = getContext();
		AssemblyEvidenceSource aes = AES(pc);
		StructuralVariationCallBuilder cb = new StructuralVariationCallBuilder(pc, (VariantContextDirectedEvidence)minimalBreakend()
				.breakpoint(BP, "GT").make());
		cb.addEvidence(new rsc(1, true));
		cb.addEvidence(new rsc(2, true));
		cb.addEvidence(new rsc(3, false));
		cb.addEvidence(new rrsc(1, false));
		cb.addEvidence(new rrsc(7, false));
		
		List<DirectedEvidence> support = Lists.<DirectedEvidence>newArrayList(
				new rsc(4, true)); 
		SAMRecord ass = AssemblyFactory.createAnchoredBreakend(pc, aes, new SequentialIdGenerator("asm"), BP.direction, support, BP.referenceIndex, BP.end, 1, B("TT"), B("TT"));
		cb.addEvidence(incorporateRealignment(AES(), ass, ImmutableList.of(withMapq(44, onNegative(Read(BP.referenceIndex2, BP.start2, "1M")))[0])));
	
		support = Lists.<DirectedEvidence>newArrayList(
				new rsc(3, false),
				new rsc(7, false),  // not supporting
				new rsc(8, false));  // not supporting
		ass = AssemblyFactory.createAnchoredBreakend(pc, aes, new SequentialIdGenerator("asm"), BP.direction, support, BP.referenceIndex2, BP.end2, 1, B("TT"), B("TT"));
		SAMRecord remote = withMapq(44, onNegative(Read(BP.referenceIndex, BP.start, "1M")))[0];
		incorporateRealignment(aes, ass, ImmutableList.of(remote));
		cb.addEvidence(asAssemblyEvidence(remote));
		
		VariantContextDirectedBreakpoint bp = (VariantContextDirectedBreakpoint)cb.make();
		
		assertAttr(VcfInfoAttributes.BREAKPOINT_ASSEMBLY_READ_COUNT, VcfFormatAttributes.BREAKPOINT_ASSEMBLY_READ_COUNT, new int[] { 3, 1, }, bp);
	}
	@Test
	@Ignore("Enchancement")
	public void how_should_we_count_spanning_assemblies() {
		fail();
	}
	/**
	 * We do this to prevent inflation of support for a breakpoint by an assembly that
	 * doesn't not actually support that breakend 
	 */
	@Test
	public void should_exclude_breakpoint_evidence_contributing_to_breakend_assembly_but_not_supporting_breakpoint() {
		StructuralVariationCallBuilder cb = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)minimalBreakend()
				.breakpoint(BP, "GT").make());
		List<DirectedEvidence> support = Lists.<DirectedEvidence>newArrayList(new rsc_not_supporting_breakpoint(1, true)); 
		SAMRecord ass = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), BP.direction, support, BP.referenceIndex, BP.end, 1, B("TT"), B("TT"));
		cb.addEvidence(asAssemblyEvidence(ass));
		VariantContextDirectedBreakpoint call = (VariantContextDirectedBreakpoint)cb.make();
		Assert.assertArrayEquals(new int[] { 0,  0, }, asIntTN(call.getAttribute(VcfInfoAttributes.BREAKPOINT_SPLITREAD_COUNT.attribute())));
		Assert.assertArrayEquals(new int[] { 0,  0, }, asIntTN(call.getAttribute(VcfInfoAttributes.BREAKEND_SOFTCLIP_COUNT.attribute())));
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
		Assert.assertEquals(1, (int)v.getAttribute(VcfInfoAttributes.BREAKEND_SOFTCLIP_COUNT.attribute()));
		Assert.assertEquals(3, (int)v.getAttribute(VcfInfoAttributes.BREAKPOINT_SPLITREAD_COUNT.attribute()));
		Assert.assertEquals(1, (int)v.getAttribute(VcfInfoAttributes.BREAKEND_UNMAPPEDMATE_COUNT.attribute()));
		Assert.assertEquals(1, (int)v.getAttribute(VcfInfoAttributes.BREAKPOINT_READPAIR_COUNT.attribute()));
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
		tes.category = 1;
		StructuralVariationCallBuilder cb = new StructuralVariationCallBuilder(pc, (VariantContextDirectedEvidence)minimalBreakend()
				.breakend(new BreakendSummary(0, BWD, 1, 1, 10), "").make());
		cb.addEvidence(SCE(BWD, nes, Read(0, 5, "1S2M")));
		cb.addEvidence(SCE(BWD, tes, Read(0, 4, "3S4M")));
		cb.addEvidence(SCE(BWD, tes, Read(0, 3, "5S6M")));
		cb.addEvidence(NRRP(nes, DP(0, 11, "5M", false, 1, 100, "6M", true)));
		cb.addEvidence(NRRP(nes, DP(0, 12, "6M", false, 1, 100, "7M", true)));
		cb.addEvidence(NRRP(tes, OEA(0, 13, "8M", false)));
		cb.addEvidence(NRRP(nes, OEA(0, 15, "9M", false)));
		cb.addEvidence(ass1());
		cb.addEvidence(ass2());
		cb.referenceReads(new int[] { 7, 8} );
		cb.referenceSpanningPairs(new int[] { 9, 10 } );
		return cb.make();
	}
	public SoftClipEvidence ass1() {
		ProcessingContext pc = getContext();
		MockSAMEvidenceSource nes = new MockSAMEvidenceSource(pc);
		MockSAMEvidenceSource tes = new MockSAMEvidenceSource(pc);
		tes.category = 1;
		List<DirectedEvidence> support = Lists.newArrayList();
		support.add(SCE(BWD, nes, Read(0, 10, "4S1M")));
		support.add(SCE(BWD, tes, Read(0, 10, "3S5M")));
		support.add(SCE(BWD, tes, Read(0, 10, "3S6M")));
		support.add(NRRP(nes, OEA(0, 15, "5M", false)));
		support.add(NRRP(tes, OEA(0, 16, "5M", false)));
		support.add(NRRP(tes, OEA(0, 17, "5M", false)));
		support.add(NRRP(tes, DP(0, 1, "2M", true, 0, 15, "5M", false)));
		support.add(NRRP(tes, DP(0, 2, "2M", true, 0, 16, "5M", false)));
		support.add(NRRP(nes, DP(0, 3, "2M", true, 0, 17, "10M", false)));
		SAMRecord ass = AssemblyFactory.createAnchoredBreakend(pc, AES(), new SequentialIdGenerator("asm"), BWD, support, 0, 10, 5, B("CGTAAAAT"), new byte[] { 0,1,2,3,4,5,6,7});
		return (SoftClipEvidence)asAssemblyEvidence(ass);
	}
	public SplitReadEvidence ass2() {
		ProcessingContext pc = getContext();
		MockSAMEvidenceSource nes = new MockSAMEvidenceSource(pc);
		MockSAMEvidenceSource tes = new MockSAMEvidenceSource(pc);
		tes.category = 1;
		List<DirectedEvidence> support = Lists.newArrayList();
		support.add(SCE(BWD, nes, Read(0, 10, "4S1M")));
		support.add(SCE(BWD, tes, Read(0, 10, "3S5M")));
		support.add(SCE(BWD, tes, Read(0, 10, "3S6M")));
		support.add(NRRP(nes, OEA(0, 15, "5M", false)));
		support.add(NRRP(tes, OEA(0, 16, "5M", false)));
		support.add(NRRP(tes, OEA(0, 17, "5M", false)));
		support.add(NRRP(tes, DP(0, 1, "2M", true, 0, 15, "5M", false)));
		support.add(NRRP(tes, DP(0, 2, "2M", true, 0, 16, "5M", false)));
		support.add(NRRP(nes, DP(0, 3, "2M", true, 0, 17, "10M", false)));
		SAMRecord e = AssemblyFactory.createAnchoredBreakend(pc, AES(), new SequentialIdGenerator("asm"), BWD, support, 0, 10, 5, B("GGTAAAAC"), new byte[] { 7,6,5,4,3,2,1,0});
		SAMRecord ra = Read(1, 102, "1S1M1S");
		ra.setReadBases(B("GGT"));
		ra.setMappingQuality(45);
		ra.setBaseQualities(new byte[] { 0,1,2});
		SingleReadEvidence ev = incorporateRealignment(AES(), e, ImmutableList.of(ra));
		return (SplitReadEvidence)ev;
	}
	@Test(expected=IllegalArgumentException.class)
	public void evidence_must_support_call() {
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(),
				(VariantContextDirectedEvidence)minimalBreakend()
					.breakend(new BreakendSummary(0, FWD, 1), null)
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
	@Ignore("need better model")
	public void should_somatic_if_tumour_BAF_greater_than_normal() {
		StructuralVariationCallBuilder cb = cb(
				new sc(1, true),
				new sc(2, true),
				new sc(3, true),
				new sc(4, true));
		cb.referenceReads(new int[] { 10, 2 });
		cb.referenceSpanningPairs(new int[] { 10, 2 });
		assertTrue(cb.make().hasAttribute(VCFConstants.SOMATIC_KEY));
	}
	@Test
	@Ignore("TODO: NYI: Not Yet Implemented")
	public void somatic_p_value_should_be_calculated_from_coverage_at_both_ends_of_the_breakend() {
		fail();
	}
	@Ignore("need better model")
	@Test
	public void should_call_somatic_from_assembly_evidence() {
		ProcessingContext pc = getContext();
		MockSAMEvidenceSource tes = new MockSAMEvidenceSource(pc);
		tes.category = 1;
		List<DirectedEvidence> support = Lists.newArrayList();
		support.add(SCE(BWD, tes, Read(0, 10, "3S5M")));
		support.add(SCE(BWD, tes, Read(0, 10, "3S6M")));
		support.add(NRRP(tes, OEA(0, 16, "5M", false)));
		support.add(NRRP(tes, OEA(0, 17, "5M", false)));
		support.add(NRRP(tes, DP(0, 10, "2M", false, 1, 100, "5M", false)));
		support.add(NRRP(tes, DP(0, 10, "2M", false, 1, 101, "5M", false)));
		SAMRecord ass = AssemblyFactory.createAnchoredBreakend(pc, AES(), new SequentialIdGenerator("asm"), BWD, support, 0, 10, 5, B("CGTAAAAT"), new byte[] { 0,1,2,3,4,5,6,7});
		support.add(asAssemblyEvidence(ass));
		VariantContextDirectedEvidence e = (VariantContextDirectedEvidence)cb(support.toArray(new DirectedEvidence[0]))
			.referenceReads(new int[] { 10, 10 })
			.referenceSpanningPairs(new int[] { 10, 10 })
			.make();
		assertTrue(e.hasAttribute(VCFConstants.SOMATIC_KEY));
	}
	@Test
	public void only_read_pair_should_set_imprecise_header() {
		assertTrue(CallSV(NRRP(OEA(0, 1, "1M", true))).hasAttribute(VcfSvConstants.IMPRECISE_KEY));
		// Should still be inexact event if our CI is 1bp because we don't know if there
		// is any untemplated sequence included in the breakpoint
		assertTrue(CallSV(NRRP(SES(300, 300), DP(0, 1, "1M", true, 0, 2, "1M", false))).hasAttribute(VcfSvConstants.IMPRECISE_KEY));
	}
	@Test
	public void unanchored_assembly_should_set_imprecise_header() {
		SAMRecord ass = AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), new BreakendSummary(0, FWD, 1, 1, 2), null, B("A"), B("A"), new int[] {0, 0});
		VariantContextDirectedEvidence dba = CallSV(ass);
		assertTrue(dba.hasAttribute(VcfSvConstants.IMPRECISE_KEY));
	}
	@Test
	public void should_set_HOMLEN_HOMSEQ_for_microhomology() {
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
			breakpoint(new BreakpointSummary(0, FWD, 10, 10, 20, 0, BWD, 40, 30, 40), "");
			phredScore(20);
		}}.make());
		List<DirectedEvidence> evidence = new ArrayList<DirectedEvidence>();
		SAMRecord r = withAttr("ez", "", withAttr("SA", "polyA,40,+,5S5M,0,0", Read(0, 11, "5M5S")))[0];
		evidence.add(asAssemblyEvidence(r));
		for (DirectedEvidence e : evidence) {
			builder.addEvidence(e);
		}
		VariantContextDirectedEvidence call = builder.make();
		assertEquals("AAAAAAAAAA", call.getAttribute(VcfSvConstants.HOMOLOGY_SEQUENCE_KEY));
		assertEquals(10, call.getAttribute(VcfSvConstants.HOMOLOGY_LENGTH_KEY));
	}
	@Test
	@Ignore() //  TODO: enchancement
	public void should_use_INS_UNKNOWN_symbolic_allele_for_unknown_untemplated_sequence() {
		// <INS:UNKNOWN> insertion of unknown sequence and length
		// CIINSLEN confidence interval around insertion  
		Assert.fail();
	}
	@Test
	@Ignore() //  TODO: enchancement
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
			breakpoint(new BreakpointSummary(0, FWD, 14, 12, 14, 1, BWD, 11, 11, 11), "GTAC");
			phredScore(20);
		}}.make());
		List<DirectedEvidence> evidence = new ArrayList<DirectedEvidence>();
		evidence.add(incorporateRealignment(AES(),
				AssemblyFactory.createAnchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), FWD, null, 0, 12, 1, B("NTTTT"), B("NTTTT")),
				ImmutableList.of(new SAMRecord(getContext().getBasicSamHeader()) {{
					setReferenceIndex(1);
					setAlignmentStart(11);
					setCigarString("4M");
					setReadNegativeStrandFlag(false);
					setReadBases(B("TTTT"));
					setMappingQuality(21);
					}})));		
		evidence.add(SR(Read(0, 12, "1M3S"), Read(1, 11, "3M")));
		for (DirectedEvidence e : evidence) {
			builder.addEvidence(e);
		}
		VariantContextDirectedEvidence call = builder.make();
		assertEquals(new BreakpointSummary(0, FWD, 12, 12, 12, 1, BWD, 11, 11, 11), call.getBreakendSummary());
	}
	@Test
	public void should_adjust_call_bounds_based_on_best_soft_clip() {
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
			breakpoint(new BreakpointSummary(0, FWD, 14, 12, 14, 1, BWD, 11, 11, 11), "GTAC");
			phredScore(20);
		}}.make());
		List<DirectedEvidence> evidence = new ArrayList<DirectedEvidence>();
		evidence.add(SR(withSequence("NNNN", Read(0, 12, "1M3S"))[0], withSequence("NNN", withMapq(14, Read(1, 11, "3M")))[0]));
		for (DirectedEvidence e : evidence) {
			builder.addEvidence(e);
		}
		VariantContextDirectedEvidence call = builder.make();
		assertEquals(new BreakpointSummary(0, FWD, 12, 12, 12, 1, BWD, 11, 11, 11), call.getBreakendSummary());
	}
	@Test
	public void should_set_CALLED_QUAL_to_parent_quality_score() {
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
			breakpoint(new BreakpointSummary(0, FWD, 14, 12, 14, 1, BWD, 11, 11, 11), "GTAC");
			phredScore(31);
		}}.make());
		List<DirectedEvidence> evidence = new ArrayList<DirectedEvidence>();
		evidence.add(SR(Read(0, 12, "1M3S"), withMapq(14, Read(1, 11, "3M"))[0]));
		for (DirectedEvidence e : evidence) {
			builder.addEvidence(e);
		}
		VariantContextDirectedEvidence call = builder.make();
		assertEquals(31, AttributeConverter.asDouble(call.getAttribute(VcfInfoAttributes.CALLED_QUAL.attribute()), 0), 0);
	}
	@Test
	public void should_set_exact_soft_clip_bounds() {
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
			breakpoint(new BreakpointSummary(0, FWD, 11, 11, 12, 0, BWD, 10, 10, 10), "");
			phredScore(10);
		}}.make());
		builder.addEvidence(SR(withSequence("NNNN", Read(0, 12, "1M3S"))[0], withSequence("NNN", Read(0, 10, "3M"))[0]));
		VariantContextDirectedEvidence de = builder.make();
		assertEquals(new BreakpointSummary(0, FWD, 12, 12, 12, 0, BWD, 10, 10, 10), de.getBreakendSummary());
	}
	@Test(expected=IllegalArgumentException.class)
	public void indels_should_have_zero_margin_applied() {
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
			breakpoint(new BreakpointSummary(0, FWD, 11, 11, 11, 0, BWD, 12, 12, 12), "");
			phredScore(10);
		}}.make());
		builder.addEvidence(IE(Read(0, 11, "1M2I1M")));
		builder.addEvidence(IE(Read(0, 13, "1M2I1M")));
		VariantContextDirectedBreakpoint de = (VariantContextDirectedBreakpoint) builder.make();
		assertEquals(1, de.getBreakpointEvidenceCount());
	}
	@Test(expected=IllegalArgumentException.class)
	public void should_exclude_unsupporting_realigned_soft_clip() {
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
			breakpoint(new BreakpointSummary(0, FWD, 11, 11, 12, 0, BWD, 10, 10, 10), "");
			phredScore(10);
		}}.make());
		builder.addEvidence(SR(withSequence("NNNN", Read(0, 12, "1M3S"))[0], withSequence("NNN", Read(1, 10, "3M"))[0]));
	}
	@Test(expected=IllegalArgumentException.class)
	public void should_exclude_unsupporting_discordant_read_pair() {
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
			breakpoint(new BreakpointSummary(0, FWD, 11, 11, 12, 0, BWD, 10, 10, 10), "");
			phredScore(10);
		}}.make());
		builder.addEvidence(NRRP(DP(0, 10, "1M", true, 0, 2, "1M", false)));
	}
	@Test
	public void anchor_cigar_should_use_X_for_exact() {
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
			breakpoint(new BreakpointSummary(0, FWD, 11, 11, 11, 0, BWD, 10, 10, 10), "");
			phredScore(10);
		}}.make());
		VariantContextDirectedEvidence vc = builder.make();
		assertEquals("1X", vc.getAttribute(VcfInfoAttributes.SUPPORT_CIGAR.attribute()));
	}
	@Test
	public void anchor_cigar_should_use_2X_for_single_bp_imprecision() {
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
			breakpoint(new BreakpointSummary(0, FWD, 11, 11, 12, 0, BWD, 10, 10, 10), "");
			phredScore(10);
		}}.make());
		VariantContextDirectedEvidence vc = builder.make();
		assertEquals("2X", vc.getAttribute(VcfInfoAttributes.SUPPORT_CIGAR.attribute()));
	}
	@Test
	public void anchor_cigar_should_use_xnx_for_large_imprecision() {
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
			breakpoint(new BreakpointSummary(0, FWD, 10, 10, 15, 0, BWD, 10, 10, 10), "");
			phredScore(10);
		}}.make());
		VariantContextDirectedEvidence vc = builder.make();
		assertEquals("1X4N1X", vc.getAttribute(VcfInfoAttributes.SUPPORT_CIGAR.attribute()));
	}
	@Test
	public void anchor_cigar_should_include_anchoring_bases_fwd() {
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
			breakpoint(new BreakpointSummary(0, FWD, 10, 10, 15, 1, FWD, 100, 100, 100), "");
			phredScore(10);
		}}.make());
		builder.addEvidence(SR(withSequence("NNNNNNNNNNNNNNN", Read(0, 5, "8M7S"))[0], onNegative(Read(1, 100-6, "7M"))[0]));
		builder.addEvidence(NRRP(DP(0, 2, "2M", true, 1, 100, "1M", true)));
		VariantContextDirectedEvidence vc = builder.make();
		// 12345678901234567890
		//     MMMMMMMMSSSSSSS
		//  MM
		//            |>>>
		//  mmdmmmmmmmX
		//assertEquals(new BreakpointSummary(0, FWD, 10, 10, 15, 1, FWD, 100, 100, 100), vc.getBreakendSummary());
		assertEquals("2M1D7M1X", vc.getAttribute(VcfInfoAttributes.SUPPORT_CIGAR.attribute()));
	}
	@Test
	public void anchor_cigar_should_include_anchoring_bases_bwd() {
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
			breakpoint(new BreakpointSummary(0, BWD, 10, 10, 15, 0, FWD, 100, 100, 100), "");
			phredScore(10);
		}}.make());
		builder.addEvidence(NRRP(DP(0, 17, "2M", false, 0, 100, "1M", true)));
		builder.addEvidence(NRRP(DP(0, 18, "2M", false, 0, 100, "1M", true)));
		builder.addEvidence(NRRP(DP(0, 19, "2M", false, 0, 100, "1M", true)));
		builder.addEvidence(NRRP(DP(0, 25, "1M", false, 0, 100, "1M", true)));
		builder.addEvidence(NRRP(DP(0, 26, "1M", false, 0, 100, "1M", true)));
		builder.addEvidence(NRRP(DP(0, 28, "1M", false, 0, 100, "1M", true)));
		VariantContextDirectedEvidence vc = builder.make();
		//          1         2         3         4         5
		// 123456789012345678901234567890123456789012345678901234567890
		//          XNNNNX MM      M  M
		//                  MM      M
		//                   MM      
		// =        XNNNNXDMMMMDDDDMMDM
		assertEquals("1X4N1X1D4M4D2M1D1M", vc.getAttribute(VcfInfoAttributes.SUPPORT_CIGAR.attribute()));
	}
	@Test
	public void anchor_cigar_should_use_local_coordinates() {
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
			breakpoint(new BreakpointSummary(0, FWD, 100, 0, BWD, 10), "");
			phredScore(10);
		}}.make());
		builder.addEvidence(SR(withSequence("AAAAAAAAANNAAAAAAAAA", Read(0, 91, "10M10S"))[0], Read(0, 10, "10M")));
		VariantContextDirectedEvidence vc = builder.make();
		assertEquals("9M1X", vc.getAttribute(VcfInfoAttributes.SUPPORT_CIGAR.attribute()));
	}
	@Test
	@Ignore("The assembler doesn't know which orientation it was assembling.")
	public void spanning_assemblies_should_use_original_parent_assembly_direction_to_determine_local_remote_status() {
		IndelEvidence r = IE(withMapq(40, Read(2, 90, "10M1I10M"))[0]);
		ImmutableList<DirectedEvidence> rid = ImmutableList.of(r);
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
			breakpoint(new BreakpointSummary(2, FWD, 100, 2, BWD, 101), "");
			phredScore(10);
		}}.make());
		String seq = S(RANDOM).substring(100-10, 100) + "N" + S(RANDOM).substring(100, 100+10);
		ProcessingContext pc = getContext();
		builder.addEvidence(asAssemblyEvidence(AssemblyFactory.createAnchoredBreakpoint(pc, AES(pc), new SequentialIdGenerator("asm"), rid, 2, 100, 10, 2, 101, 10, B(seq), B(seq))));
		
		VariantContextDirectedBreakpoint vc = (VariantContextDirectedBreakpoint) builder.make();
		assertEquals(0, vc.getBreakpointEvidenceCountLocalAssembly());
		assertEquals(1, vc.getBreakpointEvidenceCountRemoteAssembly());
	}
	@Test
	@Ignore("Issue #17: Output breakpoint assembly sequences")
	public void breakpoint_assembly_should_be_written() {
		ProcessingContext pc = getContext();
		String seq = "CATTAATCGCAAGAGCGGGTTGTATTCGcCGCCAAGTCAGCTGAAGCACCATTACCCGAtCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGATTTTGTTTACAGCCTGTCTTATATCCTGAATAACGCACCGCCTATTCG";
		int anchor = 78;
		SAMRecord ass = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), FWD,
				null,
				6, 78, anchor,
				B(seq),
				B(40,seq.length()));
		SingleReadEvidence ae = incorporateRealignment(AES(), ass, ImmutableList.of(withQual(B(40, seq.length() - anchor), Read(2, anchor + 1, String.format("%dM", seq.length() - anchor)))[0]));
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(pc, (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
			breakpoint(new BreakpointSummary(6, FWD, 78, 2, BWD, 79), "");
			phredScore(50);
		}}.make());
		builder.addEvidence(ae);
		VariantContextDirectedEvidence e = builder.make();
		FastqRecord fr = e.getBreakendAssemblyFastq();
		assertEquals(e.getEvidenceID() + "#" + Integer.toString(anchor), fr.getReadHeader());
		assertEquals(seq, fr.getReadString());
		assertEquals(fr.getReadHeader(), e.getAttribute(VcfSvConstants.BREAKPOINT_ID_KEY));
	}
	@Test
	public void breakpoint_assembly_should_be_local() {
		ProcessingContext pc = getContext();
		String seq = "CATTAATCGCAAGAGCGGGTTGTATTCGcCGCCAAGTCAGCTGAAGCACCATTACCCGAtCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGATTTTGTTTACAGCCTGTCTTATATCCTGAATAACGCACCGCCTATTCG";
		int anchor = 78;
		SAMRecord ass = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), FWD,
				null,
				6, 78, anchor,
				B(seq),
				B(40,seq.length()));
		SAMRecord remote = withQual(B(40, seq.length() - anchor), Read(2, anchor + 1, String.format("%dM", seq.length() - anchor)))[0];
		incorporateRealignment(AES(), ass, remote);
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(pc, (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
			breakpoint(new BreakpointSummary(2, BWD, 79, 6, FWD, 78), "");
			phredScore(50);
		}}.make());
		builder.addEvidence(asAssemblyEvidence(remote));
		VariantContextDirectedEvidence e = builder.make();
		FastqRecord fr = e.getBreakendAssemblyFastq();
		assertEquals(null, fr);
		assertEquals(null, e.getAttribute(VcfSvConstants.BREAKPOINT_ID_KEY));
	}
	@Test
	@Ignore("Issue #18: Output breakend support interval")
	public void support_interval_should_be_written() {
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
			breakpoint(new BreakpointSummary(0, FWD, 100, 1, BWD, 200), "");
			phredScore(50);
		}}.make());
		builder.addEvidence(NRRP(DP(0, 96, "3M", true, 1, 200, "5M", false)));
		VariantContextDirectedEvidence e = builder.make();
		//          1        2
		//          0        0
		// 1234567890        0123456789
		//      MMM  >      <MMMMM
		//
		assertEquals(-4, e.getAttributeIntOffset(VcfInfoAttributes.SUPPORT_INTERVAL, 0));
		assertEquals(-2, e.getAttributeIntOffset(VcfInfoAttributes.SUPPORT_INTERVAL, 1));
		assertEquals(0, e.getAttributeIntOffset(VcfInfoAttributes.REMOTE_SUPPORT_INTERVAL, 0));
		assertEquals(4, e.getAttributeIntOffset(VcfInfoAttributes.REMOTE_SUPPORT_INTERVAL, 1));
	}
	@Test
	public void should_write_all_genotypes() {
		Assert.fail(); // TODO: check we have all zeros
	}
	@Test
	public void should_write_qual_per_sample() {
		Assert.fail(); // TODO: check we have all zeros
	}
	@Test
	public void should_prorata_assembly_qual() {
		Assert.fail(); // TODO: check we have all zeros
	}
}