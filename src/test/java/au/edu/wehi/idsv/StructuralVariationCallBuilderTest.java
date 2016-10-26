package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.variant.vcf.VCFConstants;

import java.util.ArrayList;
import java.util.List;

import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import au.edu.wehi.idsv.vcf.VcfAttributes;
import au.edu.wehi.idsv.vcf.VcfSvConstants;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;

public class StructuralVariationCallBuilderTest extends TestHelper {
	public final static BreakpointSummary BP = new BreakpointSummary(0, BWD, 10, 10, 1, BWD, 100, 100);
	public static class sc extends SoftClipEvidence {
		protected sc(int offset, boolean tumour) {
			super(SES(tumour), BWD, withSequence("NNNNNNNNNN", Read(0, 10, "5S5M"))[0]);
			this.offset = offset;
		}
		int offset;
		boolean tumour;
		@Override public int getLocalMapq() { return 1 + offset; }
		@Override public int getLocalBaseLength() { return 2 + offset; }
		@Override public float getBreakendQual() { return 16 + offset; }
		@Override public String getEvidenceID() { return "sc" + Integer.toString(offset); }
	}
	public static class rsc extends RealignedSoftClipEvidence {
		int offset;
		protected rsc(int offset, boolean tumour) {
			super(SES(tumour), BWD, withSequence("NNNNNNNNNN", Read(0, 10, "5S5M"))[0], onNegative(withSequence("NNNNN", Read(1, 100, "5M")))[0]);
			this.offset = offset;
		}
		@Override public int getLocalMapq() { return 1 + offset; }
		@Override public int getLocalBaseLength() { return 2 + offset; }
		@Override public int getRemoteMapq() { return 6 + offset; }
		@Override public float getBreakendQual() { return 111 + offset; }
		@Override public float getBreakpointQual() { return 112 + offset; }
		@Override public String getEvidenceID() { return "rsc" + Integer.toString(offset); }
	}
	public static class rsc_not_supporting_breakpoint extends RealignedSoftClipEvidence {
		int offset;
		protected rsc_not_supporting_breakpoint(int offset, boolean tumour) {
			super(SES(tumour), BWD, Read(0, 10, "5S5M"), onNegative(Read(2, 100, "5M"))[0]);
			this.offset = offset;
		}
		@Override public int getLocalMapq() { return 1 + offset; }
		@Override public int getLocalBaseLength() { return 2 + offset; }
		@Override public int getRemoteMapq() { return 6 + offset; }
		@Override public float getBreakendQual() { return 111 + offset; }
		@Override public float getBreakpointQual() { return 112 + offset; }
		@Override public String getEvidenceID() { return "rsc" + Integer.toString(offset); }
	}
	public static class rrsc extends RealignedRemoteSoftClipEvidence {
		int offset;
		protected rrsc(int offset, boolean tumour) {
			super(new RealignedSoftClipEvidence(SES(tumour), BWD, Read(1, 100, "5S5M"), onNegative(Read(0, 10, "5M"))[0]));
			this.offset = offset;
		}
		@Override public int getLocalMapq() { return 1 + offset; }
		@Override public int getLocalBaseLength() { return 2 + offset; }
		@Override public int getRemoteMapq() { return 6 + offset; }
		@Override public float getBreakendQual() { return 211 + offset; }
		@Override public float getBreakpointQual() { return 212 + offset; }
		@Override public String getEvidenceID() { return "Rrsc" + Integer.toString(offset); }
		@Override
		public RealignedSoftClipEvidence asLocal() { return new rrsclocal(); }
		public class rrsclocal extends RealignedSoftClipEvidence {
			protected rrsclocal() {
				super(rrsc.this.getEvidenceSource(), FWD, rrsc.this.getSAMRecord(), rrsc.this.getRealignedSAMRecord());
			}
			@Override public int getLocalMapq() { return rrsc.this.getRemoteMapq(); }
			@Override public int getRemoteMapq() { return rrsc.this.getLocalMapq(); }
			@Override public float getBreakendQual() { return rrsc.this.getBreakendQual(); }
			@Override public float getBreakpointQual() { return rrsc.this.getBreakpointQual(); }
			@Override public String getEvidenceID() { return rrsc.this.getEvidenceID().substring(1); }
		}
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
		@Override public int getLocalBaseLength() { return 2 + offset; }
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
		@Override public int getLocalBaseLength() { return 2 + offset; }
		@Override public int getRemoteMapq() { return 6 + offset; }
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
		SAMRecordAssemblyEvidence ass = AssemblyFactory.createAnchoredBreakend(pc, aes, BP.direction, Lists.transform(support, EID), BP.referenceIndex, BP.end, 1, B("TT"), B("TT"));
		ass.hydrateEvidenceSet(support);
		ass.annotateAssembly();
		cb.addEvidence(AssemblyFactory.incorporateRealignment(getContext(), ass, ImmutableList.of(withMapq(44, onNegative(Read(BP.referenceIndex2, BP.start2, "1M")))[0])));
		
		support = Lists.<DirectedEvidence>newArrayList(new um(6, true)); 
		ass = AssemblyFactory.createAnchoredBreakend(pc, aes, BP.direction, Lists.transform(support, EID), BP.referenceIndex, BP.end, 1, B("TT"), B("TT"));
		ass.hydrateEvidenceSet(support);
		ass.annotateAssembly();
		cb.addEvidence(ass);
		
		support = Lists.<DirectedEvidence>newArrayList(new rsc(6, false));
		ass = AssemblyFactory.createAnchoredBreakend(pc, aes, BP.direction, Lists.transform(support, EID), BP.referenceIndex, BP.end, 1, B("TA"), B("TA"));
		ass.hydrateEvidenceSet(support);
		ass.annotateAssembly();
		cb.addEvidence(AssemblyFactory.incorporateRealignment(pc, ass, ImmutableList.of(withMapq(1, onNegative(Read(BP.referenceIndex2, BP.start2, "1M")))[0])));
		
		support = Lists.<DirectedEvidence>newArrayList(new sc(4, false));
		ass = AssemblyFactory.createAnchoredBreakend(pc, aes, BP.direction, Lists.transform(support, EID), BP.referenceIndex2, BP.end2, 1, B("TT"), B("TT"));
		ass.hydrateEvidenceSet(support);
		ass.annotateAssembly();
		cb.addEvidence(((RealignedSAMRecordAssemblyEvidence)AssemblyFactory.incorporateRealignment(pc, ass, ImmutableList.of(withMapq(44, onNegative(Read(BP.referenceIndex, BP.start, "1M")))[0]))).asRemote());
		
		return (VariantContextDirectedBreakpoint)cb.make();
	}
	@Test
	public void should_set_BREAKEND_SOFTCLIP() {
		Assert.assertArrayEquals(new int[] { 2,  1, }, (int[])complex_bp().getAttribute(VcfAttributes.BREAKEND_SOFTCLIP_COUNT.attribute()));
		Assert.assertArrayEquals(new float[] { 2*16 + 2+3,  1*16 + 1, }, (float[])complex_bp().getAttribute(VcfAttributes.BREAKEND_SOFTCLIP_QUAL.attribute()), 0);
	}
	@Test
	public void should_set_BREAKPOINT_SOFTCLIP() {
		Assert.assertArrayEquals(new int[] { 3,  2, }, (int[])complex_bp().getAttribute(VcfAttributes.BREAKPOINT_SPLITREAD_COUNT.attribute()));
		Assert.assertArrayEquals(new float[] { 3*112 + 3+4+5,  2*112 + 1+2, }, (float[])complex_bp().getAttribute(VcfAttributes.BREAKPOINT_SPLITREAD_QUAL.attribute()), 0);
	}
	@Test
	public void should_set_BREAKPOINT_SOFTCLIP_REMOTE() {
		Assert.assertArrayEquals(new int[] { 1,  3, }, (int[])complex_bp().getAttribute(VcfAttributes.BREAKPOINT_SPLITREAD_COUNT_REMOTE.attribute()));
		Assert.assertArrayEquals(new float[] { 1*212 + 4, 3*212 + 1+2+3, }, (float[])complex_bp().getAttribute(VcfAttributes.BREAKPOINT_SPLITREAD_QUAL_REMOTE.attribute()), 0);
	}
	@Test
	public void should_set_BREAKEND_READPAIR() {
		Assert.assertArrayEquals(new int[] { 4,  1, }, (int[])complex_bp().getAttribute(VcfAttributes.BREAKEND_UNMAPPEDMATE_COUNT.attribute()));
		Assert.assertArrayEquals(new float[] { 4*6 + 1+2+3+4,  1*6 + 5, }, (float[])complex_bp().getAttribute(VcfAttributes.BREAKEND_UNMAPPEDMATE_QUAL.attribute()), 0);
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
		Assert.assertEquals(5, (float)(Float)complex_bp().getAttribute(VcfAttributes.BREAKPOINT_ASSEMBLY_QUAL_REMOTE.attribute()), 0);
	}
	@Test
	public void should_set_VcfAttribute_BREAKPOINT_ASSEMBLY_READPAIR_COUNT_CONSCRIPTED() {
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
		SAMRecordAssemblyEvidence ass = AssemblyFactory.createAnchoredBreakend(pc, aes, BP.direction, Lists.transform(support, EID), BP.referenceIndex, BP.end, 1, B("TT"), B("TT"));
		ass.hydrateEvidenceSet(support);
		ass.annotateAssembly();
		cb.addEvidence(AssemblyFactory.incorporateRealignment(getContext(), ass, ImmutableList.of(withMapq(44, onNegative(Read(BP.referenceIndex2, BP.start2, "1M")))[0])));
	
		support = Lists.<DirectedEvidence>newArrayList(
				// none support since this is the remote breakend
				new dp(3, false),
				new dp(5, true),
				new dp(6, true));
		ass = AssemblyFactory.createAnchoredBreakend(pc, aes, BP.direction, Lists.transform(support, EID), BP.referenceIndex2, BP.end2, 1, B("TT"), B("TT"));
		ass.hydrateEvidenceSet(support);
		ass.annotateAssembly();
		cb.addEvidence(((RealignedSAMRecordAssemblyEvidence)AssemblyFactory.incorporateRealignment(pc, ass, ImmutableList.of(withMapq(44, onNegative(Read(BP.referenceIndex, BP.start, "1M")))[0]))).asRemote());
		
		VariantContextDirectedBreakpoint bp = (VariantContextDirectedBreakpoint)cb.make();
		
		Assert.assertArrayEquals(new int[] { 1, 4, }, (int[])bp.getAttribute(VcfAttributes.BREAKPOINT_ASSEMBLY_READPAIR_COUNT.attribute()));
		Assert.assertArrayEquals(new int[] { 1, 3, }, (int[])bp.getAttribute(VcfAttributes.BREAKPOINT_ASSEMBLY_CONSCRIPTED_READPAIR_COUNT.attribute()));
	}
	@Test
	public void should_set_VcfAttribute_BREAKPOINT_ASSEMBLY_SPLITREAD_COUNT() {
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
				new rsc(2, true).asRemote(), // not supporting
				new rsc(4, true)); 
		SAMRecordAssemblyEvidence ass = AssemblyFactory.createAnchoredBreakend(pc, aes, BP.direction, Lists.transform(support, EID), BP.referenceIndex, BP.end, 1, B("TT"), B("TT"));
		ass.hydrateEvidenceSet(support);
		ass.annotateAssembly();
		cb.addEvidence(AssemblyFactory.incorporateRealignment(getContext(), ass, ImmutableList.of(withMapq(44, onNegative(Read(BP.referenceIndex2, BP.start2, "1M")))[0])));
	
		support = Lists.<DirectedEvidence>newArrayList(
				new rsc(3, false).asRemote(),
				new rsc(5, true).asRemote(),
				new rsc(6, true).asRemote(),
				new rsc(7, false),  // not supporting
				new rsc(8, false));  // not supporting
		ass = AssemblyFactory.createAnchoredBreakend(pc, aes, BP.direction, Lists.transform(support, EID), BP.referenceIndex2, BP.end2, 1, B("TT"), B("TT"));
		ass.hydrateEvidenceSet(support);
		ass.annotateAssembly();
		cb.addEvidence(((RealignedSAMRecordAssemblyEvidence)AssemblyFactory.incorporateRealignment(pc, ass, ImmutableList.of(withMapq(44, onNegative(Read(BP.referenceIndex, BP.start, "1M")))[0]))).asRemote());
		
		VariantContextDirectedBreakpoint bp = (VariantContextDirectedBreakpoint)cb.make();
		
		Assert.assertArrayEquals(new int[] { 3, 4, }, (int[])bp.getAttribute(VcfAttributes.BREAKPOINT_ASSEMBLY_SPLITREAD_COUNT.attribute()));
		Assert.assertArrayEquals(new int[] { 2, 1, }, (int[])bp.getAttribute(VcfAttributes.BREAKPOINT_ASSEMBLY_CONSCRIPTED_SPLITREAD_COUNT.attribute()));
	}
	@Test
	@Ignore() //  TODO: enchancement
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
		SAMRecordAssemblyEvidence ass = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), BP.direction, Lists.transform(support, EID), BP.referenceIndex, BP.end, 1, B("TT"), B("TT"));
		ass.hydrateEvidenceSet(support);
		ass.annotateAssembly();
		cb.addEvidence(ass);
		VariantContextDirectedBreakpoint call = (VariantContextDirectedBreakpoint)cb.make();
		Assert.assertArrayEquals(new int[] { 0,  0, }, asIntTN(call.getAttribute(VcfAttributes.BREAKPOINT_SPLITREAD_COUNT.attribute())));
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
		Assert.assertArrayEquals(new int[] { 0,  1, }, (int[])v.getAttribute(VcfAttributes.BREAKPOINT_SPLITREAD_COUNT.attribute()));
		Assert.assertArrayEquals(new int[] { 0,  1, }, (int[])v.getAttribute(VcfAttributes.BREAKPOINT_SPLITREAD_COUNT_REMOTE.attribute()));
		Assert.assertArrayEquals(new int[] { 0,  1, }, (int[])v.getAttribute(VcfAttributes.BREAKEND_UNMAPPEDMATE_COUNT.attribute()));
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
		tes.category = 1;
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
		cb.referenceReads(new int[] { 7, 8} );
		cb.referenceSpanningPairs(new int[] { 9, 10 } );
		return cb.make();
	}
	public SAMRecordAssemblyEvidence ass1() {
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
		SAMRecordAssemblyEvidence ass = AssemblyFactory.createAnchoredBreakend(pc, AES(), BWD, Lists.transform(support, EID), 0, 10, 5, B("CGTAAAAT"), new byte[] { 0,1,2,3,4,5,6,7});
		ass.hydrateEvidenceSet(support);
		ass.annotateAssembly();
		return ass;
	}
	public RealignedSAMRecordAssemblyEvidence ass2() {
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
		SAMRecordAssemblyEvidence e = AssemblyFactory.createAnchoredBreakend(pc, AES(), BWD, Lists.transform(support, EID), 0, 10, 5, B("GGTAAAAC"), new byte[] { 7,6,5,4,3,2,1,0});
		e.hydrateEvidenceSet(support);
		e.annotateAssembly();
		SAMRecord ra = Read(1, 102, "1S1M1S");
		ra.setReadBases(B("GGT"));
		ra.setMappingQuality(45);
		ra.setBaseQualities(new byte[] { 0,1,2});
		return (RealignedSAMRecordAssemblyEvidence) AssemblyFactory.incorporateRealignment(getContext(), e, ImmutableList.of(ra));
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
		SAMRecordAssemblyEvidence ass = AssemblyFactory.createAnchoredBreakend(pc, AES(), BWD, Lists.transform(support, EID), 0, 10, 5, B("CGTAAAAT"), new byte[] { 0,1,2,3,4,5,6,7});
		ass.hydrateEvidenceSet(support);
		ass.annotateAssembly();
		support.add(ass);
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
		SAMRecordAssemblyEvidence ass = AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new BreakendSummary(0, FWD, 1, 2), null, B("A"), B("A"), new int[] {0, 0});
		ass.annotateAssembly();
		VariantContextDirectedEvidence dba = CallSV(ass);
		assertTrue(dba.hasAttribute(VcfSvConstants.IMPRECISE_KEY));
	}
	@Test
	public void should_set_HOMLEN_HOMSEQ_for_microhomology() {
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
			breakpoint(new BreakpointSummary(0, FWD, 10, 20, 0, BWD, 30, 40), "");
			phredScore(20);
		}}.make());
		List<DirectedEvidence> evidence = new ArrayList<DirectedEvidence>();
		evidence.add(SoftClipEvidence.create(SES(), FWD, Read(0, 11, "5M5S"), Read(0, 35, "5M")));
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
			breakpoint(new BreakpointSummary(0, FWD, 12, 14, 1, BWD, 11, 11), "GTAC");
			phredScore(20);
		}}.make());
		List<DirectedEvidence> evidence = new ArrayList<DirectedEvidence>();
		evidence.add(AssemblyFactory.incorporateRealignment(getContext(),
				AssemblyFactory.createAnchoredBreakend(getContext(), AES(), FWD, null, 0, 12, 1, B("NTTTT"), B("NTTTT")).annotateAssembly(),
				ImmutableList.of(new SAMRecord(getContext().getBasicSamHeader()) {{
					setReferenceIndex(1);
					setAlignmentStart(11);
					setCigarString("4M");
					setReadNegativeStrandFlag(false);
					setReadBases(B("TTTT"));
					setMappingQuality(21);
					}})));		
		evidence.add(SoftClipEvidence.create(SES(), FWD, Read(0, 12, "1M3S"), Read(1, 11, "3M")));
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
		evidence.add(SoftClipEvidence.create(SES(), FWD, Read(0, 12, "1M3S"), withSequence("NNN", withMapq(14, Read(1, 11, "3M")))[0]));
		evidence.add(SoftClipEvidence.create(SES(), FWD, Read(0, 13, "1M3S"), withSequence("NNN", withMapq(15, Read(1, 11, "3M")))[0]));
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
		evidence.add(SoftClipEvidence.create(SES(), FWD, Read(0, 12, "1M3S"), withMapq(14, Read(1, 11, "3M"))[0]));
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
		builder.addEvidence(SoftClipEvidence.create(SES(), FWD, withSequence("NNNN", Read(0, 12, "1M3S"))[0], withSequence("NNN", Read(0, 10, "3M"))[0]));
		VariantContextDirectedEvidence de = builder.make();
		assertEquals(new BreakpointSummary(0, FWD, 12, 12, 0, BWD, 10, 10), de.getBreakendSummary());
	}
	@Test(expected=IllegalArgumentException.class)
	public void indels_should_have_zero_margin_applied() {
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
			breakpoint(new BreakpointSummary(0, FWD, 11, 11, 0, BWD, 12, 12), "");
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
			breakpoint(new BreakpointSummary(0, FWD, 11, 12, 0, BWD, 10, 10), "");
			phredScore(10);
		}}.make());
		builder.addEvidence(SoftClipEvidence.create(SES(), FWD, withSequence("NNNN", Read(0, 12, "1M3S"))[0], withSequence("NNN", Read(1, 10, "3M"))[0]));
	}
	@Test(expected=IllegalArgumentException.class)
	public void should_exclude_unsupporting_discordant_read_pair() {
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
			breakpoint(new BreakpointSummary(0, FWD, 11, 12, 0, BWD, 10, 10), "");
			phredScore(10);
		}}.make());
		builder.addEvidence(NRRP(DP(0, 10, "1M", true, 0, 2, "1M", false)));
	}
	@Test
	public void anchor_cigar_should_use_X_for_exact() {
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
			breakpoint(new BreakpointSummary(0, FWD, 11, 11, 0, BWD, 10, 10), "");
			phredScore(10);
		}}.make());
		VariantContextDirectedEvidence vc = builder.make();
		assertEquals("1X", vc.getAttribute(VcfAttributes.SUPPORT_CIGAR.attribute()));
	}
	@Test
	public void anchor_cigar_should_use_2X_for_single_bp_imprecision() {
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
			breakpoint(new BreakpointSummary(0, FWD, 11, 12, 0, BWD, 10, 10), "");
			phredScore(10);
		}}.make());
		VariantContextDirectedEvidence vc = builder.make();
		assertEquals("2X", vc.getAttribute(VcfAttributes.SUPPORT_CIGAR.attribute()));
	}
	@Test
	public void anchor_cigar_should_use_xnx_for_large_imprecision() {
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
			breakpoint(new BreakpointSummary(0, FWD, 10, 15, 0, BWD, 10, 10), "");
			phredScore(10);
		}}.make());
		VariantContextDirectedEvidence vc = builder.make();
		assertEquals("1X4N1X", vc.getAttribute(VcfAttributes.SUPPORT_CIGAR.attribute()));
	}
	@Test
	public void anchor_cigar_should_include_anchoring_bases_fwd() {
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
			breakpoint(new BreakpointSummary(0, FWD, 10, 15, 0, FWD, 100, 100), "");
			phredScore(10);
		}}.make());
		builder.addEvidence(SCE(FWD, Read(0, 5, "8M7S")));
		builder.addEvidence(NRRP(DP(0, 2, "2M", true, 0, 100, "1M", true)));
		VariantContextDirectedEvidence vc = builder.make();
		// 12345678901234567890
		//     MMMMMMMMSSSSSSS
		//  MM
		//          XNNNNX
		assertEquals("2M1D5M1X4N1X", vc.getAttribute(VcfAttributes.SUPPORT_CIGAR.attribute()));
	}
	@Test
	public void anchor_cigar_should_include_anchoring_bases_bwd() {
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
			breakpoint(new BreakpointSummary(0, BWD, 10, 15, 0, FWD, 100, 100), "");
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
		assertEquals("1X4N1X1D4M4D2M1D1M", vc.getAttribute(VcfAttributes.SUPPORT_CIGAR.attribute()));
	}
	@Test
	public void anchor_cigar_should_use_local_coordinates() {
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
			breakpoint(new BreakpointSummary(0, FWD, 100, 100, 0, BWD, 10, 10), "");
			phredScore(10);
		}}.make());
		builder.addEvidence(((RealignedSoftClipEvidence)SCE(FWD, SES(), Read(0, 91, "10M10S"), Read(0, 10, "10S10M"))));
		builder.addEvidence(((RealignedSoftClipEvidence)SCE(BWD, SES(), Read(0, 10, "10S10M"), Read(0, 91, "10M10S"))).asRemote());
		SAMRecordAssemblyEvidence ass = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), BWD, null, 0, 10, 10, B("NNNNNNNNNNNNNNNNNNNN"), B("00000000001111111111"));
		ass.annotateAssembly();
		ass = AssemblyFactory.incorporateRealignment(getContext(), ass, ImmutableList.of(Read(0, 91, "10M")));
		ass = ((RealignedSAMRecordAssemblyEvidence)ass).asRemote();
		builder.addEvidence(ass);
		VariantContextDirectedEvidence vc = builder.make();
		assertEquals("9M1X", vc.getAttribute(VcfAttributes.SUPPORT_CIGAR.attribute()));
	}
	@Test
	public void spanning_assemblies_should_use_original_parent_assembly_direction_to_determine_local_remote_status() {
		SpannedIndelEvidence r = IE(withMapq(40, Read(2, 90, "10M1I10M"))[0]);
		ImmutableList<String> rid = ImmutableList.of(r.getEvidenceID());
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
			breakpoint(new BreakpointSummary(2, FWD, 100, 100, 2, BWD, 101, 101), "");
			phredScore(10);
		}}.make());
		String seq = S(RANDOM).substring(100-10, 100) + "N" + S(RANDOM).substring(100, 100+10);
		ProcessingContext pc = getContext();
		builder.addEvidence(AssemblyFactory.createAnchoredBreakend(pc, AES(pc), BWD, rid, 2, 101, 10, B(seq), B(seq)).realign(5, 0)
				.hydrateEvidenceSet(r).annotateAssembly().getSpannedIndels().get(0));
		builder.addEvidence(AssemblyFactory.createAnchoredBreakend(pc, AES(pc), FWD, rid, 2, 100, 10, B(seq), B(seq)).realign(5, 0)
				.hydrateEvidenceSet(r).annotateAssembly().getSpannedIndels().get(0));
		builder.addEvidence(AssemblyFactory.createAnchoredBreakend(pc, AES(pc), FWD, rid, 2, 100, 10, B(seq), B(seq)).realign(5, 0)
				.hydrateEvidenceSet(r).annotateAssembly().getSpannedIndels().get(0));
		builder.addEvidence(AssemblyFactory.createAnchoredBreakpoint(pc, AES(pc), rid, 2, 100, 10, 2, 101, 10, B(seq), B(seq))
				.hydrateEvidenceSet(r).annotateAssembly().getSpannedIndels().get(0));
		
		VariantContextDirectedBreakpoint vc = (VariantContextDirectedBreakpoint) builder.make();
		assertEquals(3, vc.getBreakpointEvidenceCountLocalAssembly());
		assertEquals(2, vc.getBreakpointEvidenceCountRemoteAssembly());
	}
	@Test
	public void should_calculate_inexact_homology() {
		ProcessingContext pc = getContext();
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(pc, (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
			breakpoint(new BreakpointSummary(2, FWD, 78, 78, 6, BWD, 79, 79), "");
			phredScore(50);
		}}.make());
		builder.addEvidence(SCE(FWD, Read(2, 78, "1M1S"), Read(6, 79, "1S1M")));
		VariantContextDirectedEvidence e = builder.make();
		assertEquals(-78, ((int[])e.getAttribute(VcfAttributes.INEXACT_HOMPOS.attribute()))[0]);
		assertEquals(300, ((int[])e.getAttribute(VcfAttributes.INEXACT_HOMPOS.attribute()))[1]);
	}
	@Test
	@Ignore("Issue #17: Output breakpoint assembly sequences")
	public void breakpoint_assembly_should_be_written() {
		ProcessingContext pc = getContext();
		String seq = "CATTAATCGCAAGAGCGGGTTGTATTCGcCGCCAAGTCAGCTGAAGCACCATTACCCGAtCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGATTTTGTTTACAGCCTGTCTTATATCCTGAATAACGCACCGCCTATTCG";
		int anchor = 78;
		SAMRecordAssemblyEvidence ass = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), FWD,
				null,
				6, 78, anchor,
				B(seq),
				B(40,seq.length()));
		ass = AssemblyFactory.incorporateRealignment(getContext(), ass, ImmutableList.of(
				withReadName(BreakpointFastqEncoding.getRealignmentFastq(ass).getReadHeader(), withSequence(B(seq.substring(anchor)), 
						withQual(B(40, seq.length() - anchor), 
								Read(2, anchor + 1, String.format("%dM", seq.length() - anchor)))))[0]));
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(pc, (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
			breakpoint(new BreakpointSummary(6, FWD, 78, 78, 2, BWD, 79, 79), "");
			phredScore(50);
		}}.make());
		ass.annotateAssembly();
		builder.addEvidence(ass);
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
		SAMRecordAssemblyEvidence ass = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), FWD,
				null,
				6, 78, anchor,
				B(seq),
				B(40,seq.length()));
		ass = AssemblyFactory.incorporateRealignment(getContext(), ass, ImmutableList.of(
				withReadName(BreakpointFastqEncoding.getRealignmentFastq(ass).getReadHeader(), withSequence(B(seq.substring(anchor)), 
						withQual(B(40, seq.length() - anchor), 
								Read(2, anchor + 1, String.format("%dM", seq.length() - anchor)))))[0]));
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(pc, (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
			breakpoint(new BreakpointSummary(2, BWD, 79, 79, 6, FWD, 78, 78), "");
			phredScore(50);
		}}.make());
		ass.annotateAssembly();
		builder.addEvidence(((RealignedSAMRecordAssemblyEvidence)ass).asRemote());
		VariantContextDirectedEvidence e = builder.make();
		FastqRecord fr = e.getBreakendAssemblyFastq();
		assertEquals(null, fr);
		assertEquals(null, e.getAttribute(VcfSvConstants.BREAKPOINT_ID_KEY));
	}
	@Test
	@Ignore("Issue #18: Output breakend support interval")
	public void support_interval_should_be_written() {
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
			breakpoint(new BreakpointSummary(0, FWD, 100, 100, 1, BWD, 200, 200), "");
			phredScore(50);
		}}.make());
		builder.addEvidence(NRRP(DP(0, 96, "3M", true, 1, 200, "5M", false)));
		VariantContextDirectedEvidence e = builder.make();
		//          1        2
		//          0        0
		// 1234567890        0123456789
		//      MMM  >      <MMMMM
		//
		assertEquals(-4, e.getAttributeIntOffset(VcfAttributes.SUPPORT_INTERVAL, 0));
		assertEquals(-2, e.getAttributeIntOffset(VcfAttributes.SUPPORT_INTERVAL, 1));
		assertEquals(0, e.getAttributeIntOffset(VcfAttributes.REMOTE_SUPPORT_INTERVAL, 0));
		assertEquals(4, e.getAttributeIntOffset(VcfAttributes.REMOTE_SUPPORT_INTERVAL, 1));
	}
}