package au.edu.wehi.idsv;

import static org.junit.Assert.*;
import htsjdk.samtools.SAMRecord;
import htsjdk.variant.vcf.VCFConstants;

import java.util.Set;

import org.junit.Test;

import au.edu.wehi.idsv.vcf.VcfAttributes;
import au.edu.wehi.idsv.vcf.VcfAttributes.Subset;

import com.google.common.collect.Sets;

public class StructuralVariationCallBuilderTest extends TestHelper {
	public class sc extends SoftClipEvidence {
		protected sc(int offset, boolean tumour) {
			super(getContext(), SES(tumour), BWD, Read(0, 10, "5S5M"));
			this.offset = offset;
		}
		int offset;
		boolean tumour;
		@Override public int getLocalMapq() { return 1 + offset; }
		@Override public int getLocalBaseLength() { return 2 + offset; }
		@Override public int getLocalBaseCount() { return 3 + offset; }
		@Override public int getLocalMaxBaseQual() { return 4 + offset; }
		@Override public int getLocalTotalBaseQual() { return 5 + offset; }
	}
	public class rsc extends RealignedSoftClipEvidence {
		int offset;
		protected rsc(int offset, boolean tumour) {
			super(getContext(), SES(tumour), BWD, Read(0, 10, "5S5M"), Read(0, 1, "5M"));
			this.offset = offset;
		}
		@Override public int getLocalMapq() { return 1 + offset; }
		@Override public int getLocalBaseLength() { return 2 + offset; }
		@Override public int getLocalBaseCount() { return 3 + offset; }
		@Override public int getLocalMaxBaseQual() { return 4 + offset; }
		@Override public int getLocalTotalBaseQual() { return 5 + offset; }
		@Override public int getRemoteMapq() { return 6 + offset; }
		@Override public int getRemoteBaseLength() { return 7 + offset; }
		@Override public int getRemoteBaseCount() { return 8 + offset; }
		@Override public int getRemoteMaxBaseQual() { return 9 + offset; }
		@Override public int getRemoteTotalBaseQual() { return 10 + offset; }
	}
	public class um extends UnmappedMateReadPair {
		int offset;
		protected um(int offset, boolean tumour) {
			super(OEA(0, 15, "5M", false)[0],
				OEA(0, 15, "5M", false)[1],
				SES(tumour));
			this.offset = offset;
		}
		@Override public int getLocalMapq() { return 1 + offset; }
		@Override public int getLocalBaseLength() { return 2 + offset; }
		@Override public int getLocalBaseCount() { return 3 + offset; }
		@Override public int getLocalMaxBaseQual() { return 4 + offset; }
		@Override public int getLocalTotalBaseQual() { return 5 + offset; }
	}
	public class dp extends DiscordantReadPair {
		int offset;
		protected dp(int offset, boolean tumour) {
			super(DP(0, 12, "6M", false, 1, 10, "7M", false)[0],
					DP(0, 12, "6M", false, 1, 10, "7M", false)[1],
					SES(tumour));
				this.offset = offset;
			}
		@Override public int getLocalMapq() { return 1 + offset; }
		@Override public int getLocalBaseLength() { return 2 + offset; }
		@Override public int getLocalBaseCount() { return 3 + offset; }
		@Override public int getLocalMaxBaseQual() { return 4 + offset; }
		@Override public int getLocalTotalBaseQual() { return 5 + offset; }
		@Override public int getRemoteMapq() { return 6 + offset; }
		@Override public int getRemoteBaseLength() { return 7 + offset; }
		@Override public int getRemoteBaseCount() { return 8 + offset; }
		@Override public int getRemoteMaxBaseQual() { return 9 + offset; }
		@Override public int getRemoteTotalBaseQual() { return 10 + offset; }
	}
	@Test(expected=IllegalArgumentException.class)
	public void should_not_allow_unsupporting_evidence() {
		StructuralVariationCallBuilder cb = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)minimalBreakend()
				.breakend(new BreakendSummary(0, BWD, 1, 10), null).make());
		cb.addEvidence(SCE(BWD, Read(1, 5, "1S2M")));
	}
	@Test
	public void should_set_VcfAttribute_LOG_LIKELIHOOD_RATIO() {
		assertTrue(big().hasAttribute(VcfAttributes.LOG_LIKELIHOOD_RATIO.attribute()));
		assertNotNull(big().getAttribute(VcfAttributes.LOG_LIKELIHOOD_RATIO.attribute()));
		assertTrue((double)(Double)big().getAttribute(VcfAttributes.LOG_LIKELIHOOD_RATIO.attribute()) > 0d);
	}
	@Test
	public void should_not_set_VcfAttribute_LOG_LIKELIHOOD_RATIO_BREAKPOINT() {
		assertFalse(big().hasAttribute(VcfAttributes.LOG_LIKELIHOOD_RATIO_BREAKPOINT.attribute()));
	}
	@Test
	public void should_set_VcfAttribute_REFERENCE_COUNT_READ() {
		assertEquals(7, big().getReferenceReadCount(Subset.NORMAL));
		assertEquals(8, big().getReferenceReadCount(Subset.TUMOUR));
	}
	@Test
	public void should_set_VcfAttribute_REFERENCE_COUNT_READPAIR() {
		assertEquals(9, big().getReferenceReadPairCount(Subset.NORMAL));
		assertEquals(10, big().getReferenceReadPairCount(Subset.TUMOUR));
	}
	@Test
	public void should_set_VcfAttribute_READPAIR_EVIDENCE_COUNT() {
		VariantContextDirectedEvidence e = b(
			new dp(0, true),
			new dp(1, true),
			new dp(2, false),
			new dp(3, false),
			new um(4, true),
			new um(5, true),
			new um(6, true),
			new um(7, false));
		assertEquals(3, e.getEvidenceCountReadPair(Subset.NORMAL));
		assertEquals(5, e.getEvidenceCountReadPair(Subset.TUMOUR));
	}
	@Test
	public void should_set_VcfAttribute_READPAIR_LOG_LIKELIHOOD_RATIO() {
		VariantContextDirectedEvidence e = b(
				new dp(0, true),
				new dp(1, true),
				new dp(2, false),
				new dp(3, false),
				new um(4, true),
				new um(5, true),
				new um(6, true),
				new um(7, false));
		assertTrue(e.getBreakendLogLikelihoodReadPair(Subset.NORMAL) > 0d);
		assertTrue(e.getBreakendLogLikelihoodReadPair(Subset.TUMOUR) > 0d);
	}
	@Test
	public void should_set_VcfAttribute_READPAIR_MAPPED_READPAIR() {
		VariantContextDirectedEvidence e = b(
				new dp(0, true),
				new dp(1, true),
				new dp(2, false),
				new um(4, true),
				new um(5, true),
				new um(6, true),
				new um(7, false));
		assertEquals(1, e.getMappedEvidenceCountReadPair(Subset.NORMAL));
		assertEquals(2, e.getMappedEvidenceCountReadPair(Subset.TUMOUR));
	}
	@Test
	public void should_set_VcfAttribute_READPAIR_MAPQ_LOCAL_MAX() {
		VariantContextDirectedEvidence e = b(
				new dp(0, true),
				new dp(1, true),
				new dp(2, false),
				new dp(3, false),
				new um(4, true),
				new um(5, true),
				new um(6, true),
				new um(7, false));
		assertEquals(7+1, e.getMapqReadPairLocalMax(Subset.NORMAL));
		assertEquals(6+1, e.getMapqReadPairLocalMax(Subset.TUMOUR));
	}
	@Test
	public void should_set_VcfAttribute_READPAIR_MAPQ_LOCAL_TOTAL() {
		VariantContextDirectedEvidence e = b(
				new dp(0, true),
				new dp(1, true),
				new dp(2, false),
				new dp(3, false),
				new um(4, true),
				new um(5, true),
				new um(6, true),
				new um(7, false));
		assertEquals(3*1 + (2+3+7), e.getMapqReadPairLocalTotal(Subset.NORMAL));
		assertEquals(5*1 + (0+1+4+5+6), e.getMapqReadPairLocalTotal(Subset.TUMOUR));
	}
	@Test
	public void should_set_VcfAttribute_READPAIR_MAPQ_REMOTE_MAX() {
		VariantContextDirectedEvidence e = b(
				new dp(0, true),
				new dp(1, true),
				new dp(2, false),
				new dp(3, false),
				new um(4, true),
				new um(5, true),
				new um(6, true),
				new um(7, false));
		assertEquals(3+6, e.getMapqReadPairRemoteMax(Subset.NORMAL));
		assertEquals(1+6, e.getMapqReadPairRemoteMax(Subset.TUMOUR));
	}
	@Test
	public void should_set_VcfAttribute_READPAIR_MAPQ_REMOTE_TOTAL() {
		VariantContextDirectedEvidence e = b(
				new dp(0, true),
				new dp(1, true),
				new dp(2, false),
				new dp(3, false),
				new um(4, true),
				new um(5, true),
				new um(6, true),
				new um(7, false));
		assertEquals(2*6 + (2+3), e.getMapqReadPairRemoteTotal(Subset.NORMAL));
		assertEquals(2*6 + (0+1), e.getMapqReadPairRemoteTotal(Subset.TUMOUR));
	}
	@Test
	public void should_set_VcfAttribute_SOFTCLIP_EVIDENCE_COUNT() {
		VariantContextDirectedEvidence e = b(
				new sc(0, true),
				new sc(1, true),
				new sc(2, false),
				new sc(3, false),
				new sc(4, false),
				new rsc(5, true),
				new rsc(6, true),
				new rsc(7, true),
				new rsc(8, false));
		assertEquals(4, e.getEvidenceCountSoftClip(Subset.NORMAL));
		assertEquals(5, e.getEvidenceCountSoftClip(Subset.TUMOUR));
	}
	@Test
	public void should_set_VcfAttribute_SOFTCLIP_LOG_LIKELIHOOD_RATIO() {
		VariantContextDirectedEvidence e = b(
				new sc(0, true),
				new sc(1, true),
				new sc(2, false),
				new sc(3, false),
				new sc(4, false),
				new rsc(5, true),
				new rsc(6, true),
				new rsc(7, true),
				new rsc(8, false));
		assertTrue(e.getBreakendLogLikelihoodSoftClip(Subset.NORMAL) > 0d);
		assertTrue(e.getBreakendLogLikelihoodSoftClip(Subset.TUMOUR) > 0d);
	}
	@Test
	public void should_set_VcfAttribute_SOFTCLIP_MAPPED() {
		VariantContextDirectedEvidence e = b(
				new sc(0, true),
				new sc(1, true),
				new sc(2, false),
				new sc(3, false),
				new sc(4, false),
				new rsc(5, true),
				new rsc(6, true),
				new rsc(7, true),
				new rsc(8, false));
		assertEquals(1, e.getMappedEvidenceCountSoftClip(Subset.NORMAL));
		assertEquals(3, e.getMappedEvidenceCountSoftClip(Subset.TUMOUR));
	}
	@Test
	public void should_set_VcfAttribute_SOFTCLIP_MAPQ_REMOTE_TOTAL() {
		VariantContextDirectedEvidence e = b(
				new sc(0, true),
				new sc(1, true),
				new sc(2, false),
				new sc(3, false),
				new sc(4, false),
				new rsc(5, true),
				new rsc(6, true),
				new rsc(7, true),
				new rsc(8, false));
		assertEquals(6*1 + 8, e.getMapqSoftClipRemoteTotal(Subset.NORMAL));
		assertEquals(6*3 + (5+6+7), e.getMapqSoftClipRemoteTotal(Subset.TUMOUR));
	}
	@Test
	public void should_set_VcfAttribute_SOFTCLIP_MAPQ_REMOTE_MAX() {
		VariantContextDirectedEvidence e = b(
				new sc(0, true),
				new sc(1, true),
				new sc(2, false),
				new sc(3, false),
				new sc(4, false),
				new rsc(5, true),
				new rsc(6, true),
				new rsc(7, true),
				new rsc(8, false));
		assertEquals(6+8, e.getMapqSoftClipRemoteMax(Subset.NORMAL));
		assertEquals(6+7, e.getMapqSoftClipRemoteMax(Subset.TUMOUR));
	}
	@Test
	public void should_set_VcfAttribute_SOFTCLIP_LENGTH_REMOTE_TOTAL() {
		assertEquals(1, big().getLengthSoftClipTotal(Subset.NORMAL));
		assertEquals(3+5, big().getLengthSoftClipTotal(Subset.TUMOUR));
	}
	@Test
	public void should_set_VcfAttribute_SOFTCLIP_LENGTH_REMOTE_MAX() {
		assertEquals(1, big().getLengthSoftClipMax(Subset.NORMAL));
		assertEquals(5, big().getLengthSoftClipMax(Subset.TUMOUR));
	}
	private VariantContextDirectedEvidence fakeass(int offset) {
		IdsvVariantContextBuilder x = minimalBreakend()
				.breakend(new BreakendSummary(0, BWD, 1, 10), null);
		x.attribute(VcfAttributes.ASSEMBLY_EVIDENCE_COUNT.attribute(), offset + 1);
		x.attribute(VcfAttributes.ASSEMBLY_LOG_LIKELIHOOD_RATIO.attribute(), (double)(offset + 2));
		x.attribute(VcfAttributes.ASSEMBLY_MAPPED.attribute(), offset + 3);
		x.attribute(VcfAttributes.ASSEMBLY_MAPQ_REMOTE_TOTAL.attribute(), offset + 4);
		
		x.attribute(VcfAttributes.ASSEMBLY_BASE_COUNT.attribute(), new int[] { offset + 6, offset + 16});
		x.attribute(VcfAttributes.ASSEMBLY_READPAIR_COUNT.attribute(), new int[] { offset + 7, offset + 17});
		x.attribute(VcfAttributes.ASSEMBLY_SOFTCLIP_COUNT.attribute(), new int[] { offset + 8, offset + 18});
		x.attribute(VcfAttributes.ASSEMBLY_SOFTCLIP_CLIPLENGTH_TOTAL.attribute(), new int[] { offset + 9, offset + 19});
		x.attribute(VcfAttributes.ASSEMBLY_MAPQ_REMOTE_MAX.attribute(), offset + 11);
		x.attribute(VcfAttributes.ASSEMBLY_LENGTH_LOCAL_MAX.attribute(), offset + 12);
		x.attribute(VcfAttributes.ASSEMBLY_LENGTH_REMOTE_MAX.attribute(), offset + 13);
		x.attribute(VcfAttributes.ASSEMBLY_READPAIR_LENGTH_MAX.attribute(), new int[] { offset + 14, offset + 24});
		x.attribute(VcfAttributes.ASSEMBLY_SOFTCLIP_CLIPLENGTH_MAX.attribute(), new int[] { offset + 15, offset + 25});
		return (VariantContextDirectedEvidence)x.make();
	}
	@Test
	public void should_sum_assembly_VcfAttributes_() {
		VariantContextDirectedEvidence e = b(fakeass(0), fakeass(100));
		assertEquals(1 + 101, e.getEvidenceCountAssembly());
		assertEquals((double)2 + 102, e.getBreakendLogLikelihoodAssembly(), 0);
		assertEquals(3 + 103, e.getMappedEvidenceCountAssembly());
		assertEquals(4 + 104, e.getMapqAssemblyRemoteTotal());
		
		assertEquals(6 + 106, e.getAssemblyBaseCount(Subset.NORMAL));
		assertEquals(16 + 116, e.getAssemblyBaseCount(Subset.TUMOUR));
		assertEquals(7 + 107, e.getAssemblySupportCountReadPair(Subset.NORMAL));
		assertEquals(17 + 117, e.getAssemblySupportCountReadPair(Subset.TUMOUR));
		assertEquals(8 + 108, e.getAssemblySupportCountSoftClip(Subset.NORMAL));
		assertEquals(18 + 118, e.getAssemblySupportCountSoftClip(Subset.TUMOUR));
		assertEquals(9 + 109, e.getAssemblySoftClipLengthTotal(Subset.NORMAL));
		assertEquals(19 + 119, e.getAssemblySoftClipLengthTotal(Subset.TUMOUR));
	}
	@Test
	public void should_max_assembly_VcfAttributes_() {
		VariantContextDirectedEvidence e = b(fakeass(0), fakeass(100));
		assertEquals(111, e.getMapqAssemblyRemoteMax());
		assertEquals(112, e.getAssemblyAnchorLengthMax());
		assertEquals(113, e.getAssemblyBreakendLengthMax());
		assertEquals(114, e.getAssemblyReadPairLengthMax(Subset.NORMAL));
		assertEquals(124, e.getAssemblyReadPairLengthMax(Subset.TUMOUR));
		assertEquals(115, e.getAssemblySoftClipLengthMax(Subset.NORMAL));
		assertEquals(125, e.getAssemblySoftClipLengthMax(Subset.TUMOUR));
	}
	@Test
	public void should_concat_assembly_strings() {
		assertEquals(2, big().getAttributeAsStringList(VcfAttributes.ASSEMBLY_CONSENSUS.attribute()).size());
		assertEquals(2, big().getAttributeAsStringList(VcfAttributes.ASSEMBLY_PROGRAM.attribute()).size());
	}
	@Test
	public void should_drop_base_quals() {
		assertFalse(big().hasAttribute(VcfAttributes.ASSEMBLY_BREAKEND_QUALS.attribute()));
	}
	public StructuralVariationCallBuilder cb(DirectedEvidence... evidence) {
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), (VariantContextDirectedEvidence)minimalBreakend()
				.breakend(new BreakendSummary(0, BWD, 1, 10), null).make());
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
				.breakend(new BreakendSummary(0, BWD, 1, 10), null).make());
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
	public VariantContextDirectedEvidence ass1() {
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
		AssemblyBuilder sb = new AssemblyBuilder(pc, AES())
			.direction(BWD)
			.anchorLength(5)
			.referenceAnchor(0, 10)
			.assemblerName("assemblerName")
			.assemblyBases(B("CGTAAAAT"))
			.assembledBaseCount(13, 45)
			.contributingEvidence(support)
			.assemblyBaseQuality(new byte[] { 0,1,2,3,4,5,6,7});
		VariantContextDirectedEvidence v = sb.makeVariant();
		return v;
	}
	public VariantContextDirectedEvidence ass2() {
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
		AssemblyBuilder sb = new AssemblyBuilder(pc, AES())
			.direction(BWD)
			.anchorLength(5)
			.referenceAnchor(0, 10)
			.assemblerName("assemblerName")
			.assemblyBases(B("GGTAAAAC"))
			.assembledBaseCount(21, 23)
			.contributingEvidence(support)
			.assemblyBaseQuality(new byte[] { 7,6,5,4,3,2,1,0});
		VariantContextDirectedEvidence v = sb.makeVariant();
		SAMRecord ra = Read(1, 102, "1S1M1S");
		ra.setReadBases(B("GGT"));
		ra.setMappingQuality(7);
		ra.setBaseQualities(new byte[] { 0,1,2});
		v = AssemblyBuilder.incorporateRealignment(getContext(), v, ra);
		return v;
	}
	@Test(expected=IllegalArgumentException.class)
	public void evidence_must_support_call() {
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(),
				(VariantContextDirectedEvidence)minimalBreakend()
					.breakend(new BreakendSummary(0, FWD, 1, 1), null)
					.make());
		builder.addEvidence(SCE(FWD, Read(0, 2, "1M5S")));
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
		support.add(NRRP(tes, DP(0, 1, "2M", true, 0, 15, "5M", false)));
		support.add(NRRP(tes, DP(0, 2, "2M", true, 0, 16, "5M", false)));
		AssemblyBuilder sb = new AssemblyBuilder(pc, AES())
			.direction(BWD)
			.anchorLength(5)
			.referenceAnchor(0, 10)
			.assemblerName("assemblerName")
			.assemblyBases(B("CGTAAAAT"))
			.assembledBaseCount(13, 45)
			.contributingEvidence(support)
			.assemblyBaseQuality(new byte[] { 0,1,2,3,4,5,6,7});
		VariantContextDirectedEvidence v = sb.makeVariant();
		VariantContextDirectedEvidence e = cb(v)
			.referenceReads(10, 10)
			.referenceSpanningPairs(10, 10)
			.make();
		assertTrue(e.hasAttribute(VCFConstants.SOMATIC_KEY));
	}
}
