package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;
import htsjdk.samtools.SAMRecord;
import htsjdk.variant.vcf.VCFConstants;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import org.junit.Ignore;
import org.junit.Test;

import au.edu.wehi.idsv.vcf.VcfAttributes;
import au.edu.wehi.idsv.vcf.VcfFilter;
import au.edu.wehi.idsv.vcf.VcfSvConstants;

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
		@Override public int getLocalMaxBaseQual() { return 4 + offset; }
		@Override public int getLocalTotalBaseQual() { return 5 + offset; }
	}
	public RSC rsc(int offset, boolean tumour) {
		try {
			return new RSC(offset, tumour);
		} catch (CloneNotSupportedException e) {
			throw new RuntimeException(e);
		}
	}
	public class RSC extends RealignedSoftClipEvidence {
		int offset;
		protected RSC(int offset, boolean tumour) throws CloneNotSupportedException {
			super(getContext(), SES(tumour), BWD, Read(0, 10, "5S5M"), Read(0, 1, "5M"));
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
	}
	public class um extends UnmappedMateReadPair {
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
	}
	public class dp extends DiscordantReadPair {
		int offset;
		protected dp(int offset, boolean tumour) {
			super(DP(0, 12, "6M", false, 1, 10, "7M", false)[0],
					DP(0, 12, "6M", false, 1, 10, "7M", false)[1],
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
		assertEquals(7, big().getReferenceReadCount(EvidenceSubset.NORMAL));
		assertEquals(8, big().getReferenceReadCount(EvidenceSubset.TUMOUR));
	}
	@Test
	public void should_set_VcfAttribute_REFERENCE_COUNT_READPAIR() {
		assertEquals(9, big().getReferenceReadPairCount(EvidenceSubset.NORMAL));
		assertEquals(10, big().getReferenceReadPairCount(EvidenceSubset.TUMOUR));
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
		assertEquals(3, e.getEvidenceCountReadPair(EvidenceSubset.NORMAL));
		assertEquals(5, e.getEvidenceCountReadPair(EvidenceSubset.TUMOUR));
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
		assertTrue(e.getBreakendLogLikelihoodReadPair(EvidenceSubset.NORMAL) > 0d);
		assertTrue(e.getBreakendLogLikelihoodReadPair(EvidenceSubset.TUMOUR) > 0d);
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
		assertEquals(1, e.getMappedEvidenceCountReadPair(EvidenceSubset.NORMAL));
		assertEquals(2, e.getMappedEvidenceCountReadPair(EvidenceSubset.TUMOUR));
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
		assertEquals(7+1, e.getMapqReadPairLocalMax(EvidenceSubset.NORMAL));
		assertEquals(6+1, e.getMapqReadPairLocalMax(EvidenceSubset.TUMOUR));
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
		assertEquals(3*1 + (2+3+7), e.getMapqReadPairLocalTotal(EvidenceSubset.NORMAL));
		assertEquals(5*1 + (0+1+4+5+6), e.getMapqReadPairLocalTotal(EvidenceSubset.TUMOUR));
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
		assertEquals(3+6, e.getMapqReadPairRemoteMax(EvidenceSubset.NORMAL));
		assertEquals(1+6, e.getMapqReadPairRemoteMax(EvidenceSubset.TUMOUR));
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
		assertEquals(2*6 + (2+3), e.getMapqReadPairRemoteTotal(EvidenceSubset.NORMAL));
		assertEquals(2*6 + (0+1), e.getMapqReadPairRemoteTotal(EvidenceSubset.TUMOUR));
	}
	@Test
	public void should_set_VcfAttribute_SOFTCLIP_EVIDENCE_COUNT() {
		VariantContextDirectedEvidence e = b(
				new sc(0, true),
				new sc(1, true),
				new sc(2, false),
				new sc(3, false),
				new sc(4, false),
				rsc(5, true),
				rsc(6, true),
				rsc(7, true),
				rsc(8, false));
		assertEquals(4, e.getEvidenceCountSoftClip(EvidenceSubset.NORMAL));
		assertEquals(5, e.getEvidenceCountSoftClip(EvidenceSubset.TUMOUR));
	}
	@Test
	public void should_set_VcfAttribute_SOFTCLIP_LOG_LIKELIHOOD_RATIO() {
		VariantContextDirectedEvidence e = b(
				new sc(0, true),
				new sc(1, true),
				new sc(2, false),
				new sc(3, false),
				new sc(4, false),
				rsc(5, true),
				rsc(6, true),
				rsc(7, true),
				rsc(8, false));
		assertTrue(e.getBreakendLogLikelihoodSoftClip(EvidenceSubset.NORMAL) > 0d);
		assertTrue(e.getBreakendLogLikelihoodSoftClip(EvidenceSubset.TUMOUR) > 0d);
	}
	@Test
	public void should_set_VcfAttribute_SOFTCLIP_MAPPED() {
		VariantContextDirectedEvidence e = b(
				new sc(0, true),
				new sc(1, true),
				new sc(2, false),
				new sc(3, false),
				new sc(4, false),
				rsc(5, true),
				rsc(6, true),
				rsc(7, true),
				rsc(8, false));
		assertEquals(1, e.getMappedEvidenceCountSoftClip(EvidenceSubset.NORMAL));
		assertEquals(3, e.getMappedEvidenceCountSoftClip(EvidenceSubset.TUMOUR));
	}
	@Test
	public void should_set_VcfAttribute_SOFTCLIP_MAPQ_REMOTE_TOTAL() {
		VariantContextDirectedEvidence e = b(
				new sc(0, true),
				new sc(1, true),
				new sc(2, false),
				new sc(3, false),
				new sc(4, false),
				rsc(5, true),
				rsc(6, true),
				rsc(7, true),
				rsc(8, false));
		assertEquals(6*1 + 8, e.getMapqSoftClipRemoteTotal(EvidenceSubset.NORMAL));
		assertEquals(6*3 + (5+6+7), e.getMapqSoftClipRemoteTotal(EvidenceSubset.TUMOUR));
	}
	@Test
	public void should_set_VcfAttribute_SOFTCLIP_MAPQ_REMOTE_MAX() {
		VariantContextDirectedEvidence e = b(
				new sc(0, true),
				new sc(1, true),
				new sc(2, false),
				new sc(3, false),
				new sc(4, false),
				rsc(5, true),
				rsc(6, true),
				rsc(7, true),
				rsc(8, false));
		assertEquals(6+8, e.getMapqSoftClipRemoteMax(EvidenceSubset.NORMAL));
		assertEquals(6+7, e.getMapqSoftClipRemoteMax(EvidenceSubset.TUMOUR));
	}
	@Test
	public void should_set_VcfAttribute_SOFTCLIP_LENGTH_REMOTE_TOTAL() {
		assertEquals(1, big().getLengthSoftClipTotal(EvidenceSubset.NORMAL));
		assertEquals(3+5, big().getLengthSoftClipTotal(EvidenceSubset.TUMOUR));
	}
	@Test
	public void should_set_VcfAttribute_SOFTCLIP_LENGTH_REMOTE_MAX() {
		assertEquals(1, big().getLengthSoftClipMax(EvidenceSubset.NORMAL));
		assertEquals(5, big().getLengthSoftClipMax(EvidenceSubset.TUMOUR));
	}
	private SAMRecordAssemblyEvidence fakeass(int offset) {
		Map<VcfAttributes, int[]> intListAttributes = new HashMap<VcfAttributes, int[]>();
		intListAttributes.put(VcfAttributes.ASSEMBLY_BASE_COUNT, new int[] { offset + 6, offset + 16});
		intListAttributes.put(VcfAttributes.ASSEMBLY_READPAIR_COUNT, new int[] { offset + 7, offset + 17});
		intListAttributes.put(VcfAttributes.ASSEMBLY_SOFTCLIP_COUNT, new int[] { offset + 8, offset + 18});
		intListAttributes.put(VcfAttributes.ASSEMBLY_SOFTCLIP_CLIPLENGTH_TOTAL, new int[] { offset + 9, offset + 19});
		intListAttributes.put(VcfAttributes.ASSEMBLY_READPAIR_LENGTH_MAX, new int[] { offset + 14, offset + 24});
		intListAttributes.put(VcfAttributes.ASSEMBLY_SOFTCLIP_CLIPLENGTH_MAX, new int[] { offset + 15, offset + 25});
		intListAttributes.put(VcfAttributes.ASSEMBLY_MAPQ_LOCAL_MAX, new int[] { offset + 16 });
		int anchorLen = offset + 12; // ASSEMBLY_LENGTH_LOCAL_MAX
		int bpLen = offset + 13; // ASSEMBLY_LENGTH_REMOTE_MAX
		byte[] b = new byte[anchorLen + bpLen];
		byte[] q = new byte[anchorLen + bpLen];
		SAMRecordAssemblyEvidence e = new SAMRecordAssemblyEvidence(null, getContext().getBasicSamHeader(), new BreakendSummary(0, BWD, 10, 10), AES(),
				anchorLen, b, q, intListAttributes);
		e.getSAMRecord().setReadName("fakeass"+Integer.toString(offset));
		SAMRecord r = new SAMRecord(getContext().getBasicSamHeader());
		r.setReferenceIndex(1);
		r.setAlignmentStart(1);
		r.setReadUnmappedFlag(false);
		r.setMappingQuality(offset + 4); // ASSEMBLY_MAPQ
		e = new RealignedSAMRecordAssemblyEvidence(getContext(), AES(), e.getSAMRecord(), r);
		return e;
	}
	@Test
	public void should_sum_assembly_attributes() {
		VariantContextDirectedEvidence e = b(fakeass(0), fakeass(100));
		assertEquals(2, e.getEvidenceCountAssembly());
		assertEquals(2, e.getMappedEvidenceCountAssembly());
		assertEquals(4 + 104, e.getMapqAssemblyTotal());
		assertEquals(6 + 106, e.getAssemblyBaseCount(EvidenceSubset.NORMAL));
		assertEquals(16 + 116, e.getAssemblyBaseCount(EvidenceSubset.TUMOUR));
		assertEquals(7 + 107, e.getAssemblySupportCountReadPair(EvidenceSubset.NORMAL));
		assertEquals(17 + 117, e.getAssemblySupportCountReadPair(EvidenceSubset.TUMOUR));
		assertEquals(8 + 108, e.getAssemblySupportCountSoftClip(EvidenceSubset.NORMAL));
		assertEquals(18 + 118, e.getAssemblySupportCountSoftClip(EvidenceSubset.TUMOUR));
		assertEquals(9 + 109, e.getAssemblySoftClipLengthTotal(EvidenceSubset.NORMAL));
		assertEquals(19 + 119, e.getAssemblySoftClipLengthTotal(EvidenceSubset.TUMOUR));
	}
	@Test
	public void should_max_assembly_attributes() {
		VariantContextDirectedEvidence e = b(fakeass(0), fakeass(100));
		assertEquals(104, e.getMapqAssemblyMax());
		assertEquals(112, e.getAssemblyAnchorLengthMax());
		assertEquals(113, e.getAssemblyBreakendLengthMax());
		assertEquals(114, e.getAssemblyReadPairLengthMax(EvidenceSubset.NORMAL));
		assertEquals(124, e.getAssemblyReadPairLengthMax(EvidenceSubset.TUMOUR));
		assertEquals(115, e.getAssemblySoftClipLengthMax(EvidenceSubset.NORMAL));
		assertEquals(125, e.getAssemblySoftClipLengthMax(EvidenceSubset.TUMOUR));
		assertEquals(116, e.getMapqAssemblyEvidenceMax());
	}
	@Test
	public void qual_should_be_llr() {
		assertEquals(big().getPhredScaledQual(), (double)AttributeConverter.asDoubleList(big().getAttribute(VcfAttributes.LOG_LIKELIHOOD_RATIO.attribute())).get(0), 0.001);
	}
	@Test
	public void should_calculate_component_llrs() {
		assertNotEquals(0d, AttributeConverter.asDoubleList(big().getAttribute(VcfAttributes.LOG_LIKELIHOOD_RATIO.attribute())).get(0));
		assertNotEquals(0d, AttributeConverter.asDoubleList(big().getAttribute(VcfAttributes.READPAIR_LOG_LIKELIHOOD_RATIO.attribute())).get(0));
		assertNotEquals(0d, AttributeConverter.asDoubleList(big().getAttribute(VcfAttributes.READPAIR_LOG_LIKELIHOOD_RATIO.attribute())).get(1));
		assertNotEquals(0d, AttributeConverter.asDoubleList(big().getAttribute(VcfAttributes.SOFTCLIP_LOG_LIKELIHOOD_RATIO.attribute())).get(0));
		assertNotEquals(0d, AttributeConverter.asDoubleList(big().getAttribute(VcfAttributes.SOFTCLIP_LOG_LIKELIHOOD_RATIO.attribute())).get(1));
		assertNotEquals(0d, AttributeConverter.asDoubleList(big().getAttribute(VcfAttributes.ASSEMBLY_LOG_LIKELIHOOD_RATIO.attribute())).get(0));
	}
	@Test
	public void should_merge_same_assembly_by_id() {
		assertEquals(1, cb(fakeass(0), fakeass(0)).make().getEvidenceCountAssembly());
	}
	@Test
	public void should_concat_assembly_strings() {
		assertEquals(2, big().getAssemblyConsensus().size());
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
		support.add(NRRP(tes, DP(0, 1, "2M", true, 0, 15, "5M", false)));
		support.add(NRRP(tes, DP(0, 2, "2M", true, 0, 16, "5M", false)));
		AssemblyEvidence ass = AssemblyFactory.createAnchored(pc, AES(), BWD,
				support, 0, 10, 5, B("CGTAAAAT"), new byte[] { 0,1,2,3,4,5,6,7}, 0, 45);
		VariantContextDirectedEvidence e = (VariantContextDirectedEvidence)cb(ass)
			.referenceReads(10, 10)
			.referenceSpanningPairs(10, 10)
			.make();
		assertTrue(e.hasAttribute(VCFConstants.SOMATIC_KEY));
	}
	@Test
	public void should_filter_small_indels() {
		assertTrue(new StructuralVariationCallBuilder(getContext(),
				(VariantContextDirectedEvidence)minimalBreakend()
					.breakpoint(new BreakpointSummary(new BreakendSummary(0, FWD, 100, 200), new BreakendSummary(0, BWD, 115, 215)), null)
					.make()).make().getFilters().contains(VcfFilter.SMALL_INDEL.filter()));
		assertFalse(new StructuralVariationCallBuilder(getContext(),
				(VariantContextDirectedEvidence)minimalBreakend()
					.breakpoint(new BreakpointSummary(new BreakendSummary(0, FWD, 100, 200), new BreakendSummary(0, BWD, 216, 316)), null)
					.make()).make().getFilters().contains(VcfFilter.SMALL_INDEL.filter()));
	}
	@Test
	public void best_assembly_should_contain_sc() {
		
	}
	@Test
	public void should_set_assembly_attribute_ASSEMBLY_EVIDENCE_COUNT() {
		assertEquals(2, big().getEvidenceCountAssembly());
	}
	@Test
	public void should_set_assembly_attribute_ASSEMBLY_LOG_LIKELIHOOD_RATIO() {
		assertTrue(big().getBreakendLogLikelihoodAssembly() > 0d);
	}
	@Test
	public void should_set_attribute_LOG_LIKELIHOOD_RATIO() {
		assertTrue(big().getBreakendLogLikelihood(null) > 0d);
	}
	@Test
	public void should_set_assembly_attribute_ASSEMBLY_MAPPED() {
		assertEquals(0, cb(new AssemblyFactoryTest().big()).make().getMappedEvidenceCountAssembly());
		assertEquals(1, cb(new AssemblyFactoryTest().bigr()).make().getMappedEvidenceCountAssembly());
		assertEquals(2, cb(new AssemblyFactoryTest().big(), fakeass(1), fakeass(2)).make().getMappedEvidenceCountAssembly());
	}
	@Test
	public void should_set_assembly_attribute_ASSEMBLY_MAPQ_REMOTE_TOTAL() {
		assertEquals(0 + 17, cb(new AssemblyFactoryTest().big(), new AssemblyFactoryTest().bigr()).make().getMapqAssemblyTotal() );
	}
	@Test
	public void mate_anchor_should_set_imprecise_header() {
		// max fragment size = 300
		VariantContextDirectedEvidence dba = CallSV(AssemblyFactory.createUnanchored(getContext(), AES(), Sets.<DirectedEvidence>newHashSet(new MockDirectedEvidence(new BreakendSummary(0, FWD, 1, 2))), B("A"), B("A"), 0, 0));
		assertTrue(dba.hasAttribute(VcfSvConstants.IMPRECISE_KEY));
	}
}
