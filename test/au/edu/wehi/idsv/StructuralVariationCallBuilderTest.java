package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import java.util.Set;

import htsjdk.samtools.SAMRecord;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;

import org.junit.Test;

import au.edu.wehi.idsv.TestHelper.MockSAMEvidenceSource;
import au.edu.wehi.idsv.vcf.VcfAttributes;

import com.google.common.collect.Sets;

public class StructuralVariationCallBuilderTest extends TestHelper {
	@Test
	public void should_group_sc_evidence_by_tumour_normal() {
		fail();
	}
	@Test
	public void should_group_rp_evidence_by_tumour_normal() {
		fail();
	}
	@Test
	public void should_group_assembly_evidence_by_tumour_normal() {
		fail();
	}
	@Test
	public void should_set_VcfAttribute_LOG_LIKELIHOOD_RATIO() {
		fail();
	}
	@Test
	public void should_set_VcfAttribute_LOG_LIKELIHOOD_RATIO_BREAKPOINT() {
		fail();
	}
	@Test
	public void should_set_VcfAttribute_REFERENCE_COUNT_READ() {
		fail();
	}
	@Test
	public void should_set_VcfAttribute_REFERENCE_COUNT_READPAIR() {
		fail();
	}
	@Test
	public void should_set_VcfAttribute_READPAIR_EVIDENCE_COUNT() {
		fail();
	}
	@Test
	public void should_set_VcfAttribute_READPAIR_LOG_LIKELIHOOD_RATIO() {
		fail();
	}
	@Test
	public void should_set_VcfAttribute_READPAIR_MAPPED_READPAIR() {
		fail();
	}
	@Test
	public void should_set_VcfAttribute_READPAIR_MAPQ_LOCAL_MAX() {
		fail();
	}
	@Test
	public void should_set_VcfAttribute_READPAIR_MAPQ_LOCAL_TOTAL() {
		fail();
	}
	@Test
	public void should_set_VcfAttribute_READPAIR_MAPQ_REMOTE_MAX() {
		fail();
	}
	@Test
	public void should_set_VcfAttribute_READPAIR_MAPQ_REMOTE_TOTAL() {
		fail();
	}
	@Test
	public void should_set_VcfAttribute_SOFTCLIP_EVIDENCE_COUNT() {
		fail();
	}
	@Test
	public void should_set_VcfAttribute_SOFTCLIP_LOG_LIKELIHOOD_RATIO() {
		fail();
	}
	@Test
	public void should_set_VcfAttribute_SOFTCLIP_MAPPED() {
		fail();
	}
	@Test
	public void should_set_VcfAttribute_SOFTCLIP_MAPQ_REMOTE_TOTAL() {
		fail();
	}
	@Test
	public void should_set_VcfAttribute_SOFTCLIP_MAPQ_REMOTE_MAX() {
		fail();
	}
	@Test
	public void should_set_VcfAttribute_SOFTCLIP_LENGTH_REMOTE_TOTAL() {
		fail();
	}
	@Test
	public void should_set_VcfAttribute_SOFTCLIP_LENGTH_REMOTE_MAX() {
		fail();
	}
	@Test
	public void should_sum_assembly_VcfAttributes_() {
		//ASSEMBLY_EVIDENCE_COUNT
		//ASSEMBLY_LOG_LIKELIHOOD_RATIO
		//ASSEMBLY_MAPPED
		//ASSEMBLY_MAPQ_REMOTE_TOTAL
		//ASSEMBLY_BASE_COUNT
		//ASSEMBLY_READPAIR_COUNT
		//ASSEMBLY_SOFTCLIP_COUNT
		//ASSEMBLY_SOFTCLIP_CLIPLENGTH_TOTAL
		fail();
	}
	@Test
	public void should_max_assembly_VcfAttributes_() {
		//ASSEMBLY_MAPQ_REMOTE_MAX
		//ASSEMBLY_LENGTH_LOCAL_MAX
		//ASSEMBLY_LENGTH_REMOTE_MAX
		//ASSEMBLY_READPAIR_LENGTH_MAX
		//ASSEMBLY_SOFTCLIP_CLIPLENGTH_MAX
		fail();
	}
	@Test
	public void should_concat_assembly_strings() {
		assertEquals(2, big().getAttributeAsStringList(VcfAttributes.ASSEMBLY_CONSENSUS.attribute()).size());
		assertEquals(2, big().getAttributeAsStringList(VcfAttributes.ASSEMBLY_PROGRAM.attribute()).size());
		assertEquals(2, big().getAttributeAsStringList(VcfAttributes.ASSEMBLY_BREAKEND_QUALS.attribute()).size());
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
		cb.addEvidence(SCE(BWD, tes, Read(0, 3, "5S46")));
		cb.addEvidence(NRRP(nes, DP(0, 11, "5M", false, 1, 100, "6M", true)));
		cb.addEvidence(NRRP(nes, DP(0, 12, "6M", false, 1, 100, "7M", true)));
		cb.addEvidence(NRRP(tes, OEA(0, 13, "8M", false)));
		cb.addEvidence(NRRP(nes, OEA(0, 15, "9M", false)));
		cb.addEvidence(ass1());
		cb.addEvidence(ass2());
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
			.referenceAnchor(0, 11)
			.assemblerName("assemblerName")
			.assemblyBases(B("CGTAAAAT"))
			.assembledBaseCount(513, 745)
			.contributingEvidence(support)
			.assemblyBaseQuality(new byte[] { 0,1,2,3,4,5,6,7});
		return sb.makeVariant();
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
			.assembledBaseCount(513, 745)
			.contributingEvidence(support)
			.assemblyBaseQuality(new byte[] { 0,1,2,3,4,5,6,7});
		VariantContextDirectedEvidence v = sb.makeVariant();
		SAMRecord ra = Read(1, 102, "1S1M1S");
		ra.setReadBases(B("GGT"));
		ra.setMappingQuality(7);
		ra.setBaseQualities(new byte[] { 0,1,2});
		return AssemblyBuilder.incorporateRealignment(getContext(), v, ra);
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
}
