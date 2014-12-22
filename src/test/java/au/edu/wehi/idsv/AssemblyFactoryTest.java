package au.edu.wehi.idsv;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;
import htsjdk.samtools.SAMRecord;

import java.util.HashSet;
import java.util.Set;

import org.junit.Test;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;


public class AssemblyFactoryTest extends TestHelper {
	@Test
	public void should_set_breakend_anchored_bwd() {
		AssemblyEvidence e = AssemblyFactory.createAnchored(
				getContext(), AES(), BWD, Sets.<DirectedEvidence>newHashSet(),
				0, 1, 2, B("GTACCC"), new byte[] { 1, 2, 3, 4, 4, 8 }, 2, 0);
		assertEquals(new BreakendSummary(0, BWD, 1, 1), e.getBreakendSummary());
	}
	@Test
	public void should_set_breakend_anchored_fwd() {
		AssemblyEvidence e = AssemblyFactory.createAnchored(
				getContext(), AES(), FWD, Sets.<DirectedEvidence>newHashSet(),
				0, 10, 2, B("GTACCC"), new byte[] { 1, 2, 3, 4, 4, 8 }, 2, 0);
		assertEquals(new BreakendSummary(0, FWD, 10, 10), e.getBreakendSummary());
	}
	@Test
	public void should_set_assembly_properties_bwd() {
		AssemblyEvidence e = AssemblyFactory.createAnchored(
				getContext(), AES(), BWD, Sets.<DirectedEvidence>newHashSet(),
				0, 1, 2, B("GTACCC"), new byte[] { 1, 2, 3, 4, 4, 8 }, 2, 0);
		assertArrayEquals(new byte[] { 1, 2, 3, 4, }, e.getBreakendQuality());
		assertEquals("GTAC", S(e.getBreakendSequence()));
		assertEquals("CC", S(e.getAssemblyAnchorSequence()));
		assertEquals("GTACCC", S(e.getAssemblySequence()));
		assertEquals(2, e.getLocalBaseLength());
		assertEquals(8, e.getLocalMaxBaseQual());
		assertEquals(4+8, e.getLocalTotalBaseQual());
	}
	@Test
	public void should_set_assembly_properties_fwd() {
		AssemblyEvidence e = AssemblyFactory.createAnchored(
				getContext(), AES(), FWD, Sets.<DirectedEvidence>newHashSet(),
				0, 1, 2, B("GTACCC"), new byte[] { 1, 3, 3, 4, 4, 8 }, 2, 0);
		assertArrayEquals(new byte[] { 3, 4, 4, 8, }, e.getBreakendQuality());
		assertEquals("ACCC", S(e.getBreakendSequence()));
		assertEquals("GT", S(e.getAssemblyAnchorSequence()));
		assertEquals("GTACCC", S(e.getAssemblySequence()));
		assertEquals(2, e.getLocalBaseLength());
		assertEquals(3, e.getLocalMaxBaseQual());
		assertEquals(1+3, e.getLocalTotalBaseQual());
	}
	@Test
	public void should_set_getAssemblySequence() {
		AssemblyEvidence e = AssemblyFactory.createAnchored(
				getContext(), AES(), FWD, Sets.<DirectedEvidence>newHashSet(),
				0, 1, 2, B("GTACCC"), new byte[] { 1, 3, 3, 4, 4, 8 }, 2, 0);
		assertEquals("GTACCC", S(e.getAssemblySequence()));
	}
	public void should_set_breakend_unanchored_fwd() {
		Set<DirectedEvidence> evidence = Sets.<DirectedEvidence>newHashSet(
				NRRP(DP(0, 10, "2M", true, 0, 10, "2M", false)));
		AssemblyEvidence e = AssemblyFactory.createUnanchored(
				getContext(), AES(), evidence,
				B("AAAAA"), B("AAAAA"), 2, 0);
		assertEquals(evidence.iterator().next().getBreakendSummary(), e.getBreakendSummary());
	}
	@Test
	public void should_set_breakend_unanchored_bwd() {
		Set<DirectedEvidence> evidence = Sets.<DirectedEvidence>newHashSet(
				NRRP(OEA(0, 1000, "5M", false)));
		AssemblyEvidence e = AssemblyFactory.createUnanchored(
				getContext(), AES(), evidence,
				B("AAAAA"), B("AAAAA"), 2, 0);
		assertEquals(evidence.iterator().next().getBreakendSummary(), e.getBreakendSummary());
	}
	@Test
	public void should_set_breakpoint_unanchored() {
		for (BreakpointSummary bs : ImmutableList.of(
				new BreakpointSummary(0, FWD, 1, 2, 3, BWD, 4, 5),
				new BreakpointSummary(0, BWD, 1, 2, 3, BWD, 4, 5)
				)) {
			AssemblyEvidence e = AssemblyFactory.createUnanchored(getContext(), AES(), Sets.<DirectedEvidence>newHashSet(new MockDirectedEvidence(bs)), B("AAAAA"), B("AAAAA"), 2, 0);
			assertEquals(bs.localBreakend(), e.getBreakendSummary());
		}
	}
	@Test
	public void should_restrict_mate_anchor_interval_based_on_anchor_positions() {
		HashSet<DirectedEvidence> evidence = Sets.<DirectedEvidence>newHashSet(
				NRRP(OEA(0, 50, "5M", true)),
				NRRP(OEA(0, 55, "5M", true)),
				NRRP(OEA(0, 60, "5M", true)),
				NRRP(OEA(0, 65, "5M", true)));
		// TODO: what to do when our assembled reads are inconsistent?
		AssemblyEvidence e = AssemblyFactory.createUnanchored(
				getContext(), AES(), evidence, B("AAAAA"), B("AAAAA"), 0, 0);
		assertEquals(Models.calculateBreakend(Lists.newArrayList(evidence)), e.getBreakendSummary());
	}
	@Test
	public void should_restrict_mate_anchor_interval_based_on_dp_interval() {
		HashSet<DirectedEvidence> evidence = Sets.<DirectedEvidence>newHashSet(
				NRRP(DP(0, 10, "1M", true, 0, 16, "1M", false)), // 10->16
				NRRP(DP(0, 12, "1M", true, 0, 19, "1M", false))); // 12->17
		AssemblyEvidence e = AssemblyFactory.createUnanchored(
				getContext(), AES(), evidence,
				B("AAAAA"), B("AAAAA"), 0, 0);
		assertEquals(new BreakendSummary(0, FWD, 12, 15), e.getBreakendSummary());
	}
	@Test
	public void should_set_assembly_attribute_ASSEMBLY_MAPQ_REMOTE_MAX() {
		assertEquals(17, bigr().getRemoteMapq());
	}
	@Test
	public void should_set_assembly_attribute_ASSEMBLY_LENGTH_LOCAL_MAX() {
		assertEquals(5, big().getAssemblyAnchorLength());
		assertEquals(5, bigr().getAssemblyAnchorLength());
	}
	@Test
	public void should_set_assembly_attribute_ASSEMBLY_LENGTH_REMOTE_MAX() {
		assertEquals(3, big().getBreakendSequence().length);
		assertEquals(3, bigr().getBreakendSequence().length);
	}
	@Test
	public void should_set_assembly_attribute_ASSEMBLY_BASE_COUNT() {	
		assertEquals(513, big().getAssemblyBaseCount(EvidenceSubset.NORMAL));
		assertEquals(745, big().getAssemblyBaseCount(EvidenceSubset.TUMOUR));
		assertEquals(513, bigr().getAssemblyBaseCount(EvidenceSubset.NORMAL));
		assertEquals(745, bigr().getAssemblyBaseCount(EvidenceSubset.TUMOUR));
	}
	@Test
	public void should_set_assembly_attribute_ASSEMBLY_READPAIR_COUNT() {
		assertEquals(2, bigr().getAssemblySupportCountReadPair(EvidenceSubset.NORMAL));
		assertEquals(4, bigr().getAssemblySupportCountReadPair(EvidenceSubset.TUMOUR));
	}
	@Test
	public void should_set_assembly_attribute_ASSEMBLY_READPAIR_LENGTH_MAX() {
		assertEquals(10, bigr().getAssemblyReadPairLengthMax(EvidenceSubset.NORMAL));
		assertEquals(5, bigr().getAssemblyReadPairLengthMax(EvidenceSubset.TUMOUR));
	}
	@Test
	public void should_set_assembly_attribute_ASSEMBLY_SOFTCLIP_COUNT() {
		assertEquals(1, bigr().getAssemblySupportCountSoftClip(EvidenceSubset.NORMAL));
		assertEquals(2, bigr().getAssemblySupportCountSoftClip(EvidenceSubset.TUMOUR));
	}
	@Test
	public void should_set_assembly_attribute_ASSEMBLY_SOFTCLIP_CLIPLENGTH_TOTAL() {	
		assertEquals(4, bigr().getAssemblySoftClipLengthTotal(EvidenceSubset.NORMAL));
		assertEquals(6, bigr().getAssemblySoftClipLengthTotal(EvidenceSubset.TUMOUR));
	}
	@Test
	public void should_set_assembly_attribute_ASSEMBLY_SOFTCLIP_CLIPLENGTH_MAX() {
		assertEquals(4, bigr().getAssemblySoftClipLengthMax(EvidenceSubset.NORMAL));
		assertEquals(3, bigr().getAssemblySoftClipLengthMax(EvidenceSubset.TUMOUR));
	}
	@Test
	public void should_set_assembly_attribute_ASSEMBLY_READPAIR_QUAL() {
		
	}
	@Test
	public void should_set_samtag_ASSEMBLY_SOFTCLIP_QUAL() {
		
	}
	@Test
	public void should_set_assembly_getLocalMapq() {
		assertEquals(10, bigr().getLocalMapq());
	}
	@Test
	public void should_set_assembly_attribute_ASSEMBLY_CONSENSUS() {
		assertEquals("CGTAAAAT", S(big().getAssemblySequence()));
		assertEquals("CGTAAAAT", S(bigr().getAssemblySequence()));
	}
	@Test
	public void should_set_assembly_attribute_ASSEMBLY_BREAKEND_QUALS() {
		assertArrayEquals(new byte[] { 0,1,2}, big().getBreakendQuality());
	}
	@Test
	public void should_set_evidence_source() {
		AssemblyEvidenceSource es = AES();
		AssemblyEvidence e = AssemblyFactory.createAnchored(
				getContext(), es, FWD, Sets.<DirectedEvidence>newHashSet(),
				0, 1, 2, B("GTACCC"), new byte[] { 1, 3, 3, 4, 4, 8 }, 2, 0);
		assertEquals(es, e.getEvidenceSource());
	}
	public SAMRecordAssemblyEvidence big() {
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
		return AssemblyFactory.createAnchored(pc, AES(), BWD, support, 0, 10, 5, B("CGTAAAAT"), new byte[] { 0,1,2,3,4,5,6,7}, 513, 745);
	}
	public RealignedSAMRecordAssemblyEvidence bigr() {
		SAMRecord ra = Read(1, 100, "1S1M1S");
		ra.setReadBases(B("CGT"));
		ra.setMappingQuality(17);
		ra.setBaseQualities(new byte[] { 0,1,2});
		return (RealignedSAMRecordAssemblyEvidence)AssemblyFactory.incorporateRealignment(getContext(), big(), ra);
	}
	@Test
	public void incorporateRealignment_should_convert_to_breakpoint() {
		assertTrue(bigr() instanceof DirectedBreakpoint);
		DirectedBreakpoint bp = (DirectedBreakpoint)bigr();
		assertEquals(1, bp.getBreakendSummary().referenceIndex2);
	}
	@Test
	public void incorporateRealignment_should_map_anchored_breakend() {
		ProcessingContext pc = getContext();
		AssemblyEvidenceSource aes = AES();
		Set<DirectedEvidence> support = Sets.newHashSet();
		
		assertEquals(new BreakpointSummary(0, FWD, 100, 100, 1, BWD, 200, 200),
				AssemblyFactory.incorporateRealignment(pc,
					AssemblyFactory.createAnchored(pc, aes, FWD, support, 0, 100, 1, B("TTT"), B("TTT"), 0, 0),
					new SAMRecord(pc.getBasicSamHeader()) {{
						setMappingQuality(40);
						setReferenceIndex(1);
						setAlignmentStart(200);
						setReadNegativeStrandFlag(false);
						setCigarString("2M");
					}}).getBreakendSummary());
		
		assertEquals(new BreakpointSummary(0, FWD, 100, 100, 1, FWD, 200, 200),
				AssemblyFactory.incorporateRealignment(pc,
					AssemblyFactory.createAnchored(pc, aes, FWD, support, 0, 100, 1, B("TTT"), B("TTT"), 0, 0),
					new SAMRecord(pc.getBasicSamHeader()) {{
						setMappingQuality(40);
						setReferenceIndex(1);
						setAlignmentStart(199);
						setReadNegativeStrandFlag(true);
						setCigarString("2M");
					}}).getBreakendSummary());
		
		assertEquals(new BreakpointSummary(0, BWD, 100, 100, 1, FWD, 200, 200),
				AssemblyFactory.incorporateRealignment(pc,
					AssemblyFactory.createAnchored(pc, aes, BWD, support, 0, 100, 1, B("TTT"), B("TTT"), 0, 0),
					new SAMRecord(pc.getBasicSamHeader()) {{
						setMappingQuality(40);
						setReferenceIndex(1);
						setAlignmentStart(199);
						setReadNegativeStrandFlag(false);
						setCigarString("2M");
					}}).getBreakendSummary());
		
		assertEquals(new BreakpointSummary(0, BWD, 100, 100, 1, BWD, 200, 200),
				AssemblyFactory.incorporateRealignment(pc,
					AssemblyFactory.createAnchored(pc, aes, BWD, support, 0, 100, 1, B("TTT"), B("TTT"), 0, 0),
					new SAMRecord(pc.getBasicSamHeader()) {{
						setMappingQuality(40);
						setReferenceIndex(1);
						setAlignmentStart(200);
						setReadNegativeStrandFlag(true);
						setCigarString("2M");
					}}).getBreakendSummary());
	}
	@Test
	public void incorporateRealignment_should_map_anchored_breakend_with_soft_clip() {
		ProcessingContext pc = getContext();
		AssemblyEvidenceSource aes = AES();
		Set<DirectedEvidence> support = Sets.newHashSet();
		
		assertEquals(new BreakpointSummary(0, FWD, 100, 100, 1, BWD, 200, 200),
				AssemblyFactory.incorporateRealignment(pc,
					AssemblyFactory.createAnchored(pc, aes, FWD, support, 0, 100, 1, B("CCCCC"), B("CCCCC"), 0, 0),
					new SAMRecord(pc.getBasicSamHeader()) {{
						setMappingQuality(40);
						setReferenceIndex(1);
						setAlignmentStart(200);
						setReadNegativeStrandFlag(false);
						setReadBases(B("CCCC"));
						setCigarString("2S2M");
					}}).getBreakendSummary());
		
		assertEquals(new BreakpointSummary(0, FWD, 100, 100, 1, FWD, 200, 200),
				AssemblyFactory.incorporateRealignment(pc,
					AssemblyFactory.createAnchored(pc, aes, FWD, support, 0, 100, 1, B("CCCCC"), B("CCCCC"), 0, 0),
					new SAMRecord(pc.getBasicSamHeader()) {{
						setMappingQuality(40);
						setReferenceIndex(1);
						setAlignmentStart(199);
						setReadNegativeStrandFlag(true);
						setReadBases(B("GGGG"));
						setCigarString("2M2S");
					}}).getBreakendSummary());
		
		assertEquals(new BreakpointSummary(0, BWD, 100, 100, 1, FWD, 200, 200),
				AssemblyFactory.incorporateRealignment(pc,
					AssemblyFactory.createAnchored(pc, aes, BWD, support, 0, 100, 1, B("CCCCC"), B("CCCCC"), 0, 0),
					new SAMRecord(pc.getBasicSamHeader()) {{
						setMappingQuality(40);
						setReferenceIndex(1);
						setAlignmentStart(199);
						setReadNegativeStrandFlag(false);
						setReadBases(B("CCCC"));
						setReadBases(B("GGGG"));
						setCigarString("2M2S");
					}}).getBreakendSummary());
		
		assertEquals(new BreakpointSummary(0, BWD, 100, 100, 1, BWD, 200, 200),
				AssemblyFactory.incorporateRealignment(pc,
					AssemblyFactory.createAnchored(pc, aes, BWD, support, 0, 100, 1, B("CCCCC"), B("CCCCC"), 0, 0),
					new SAMRecord(pc.getBasicSamHeader()) {{
						setMappingQuality(40);
						setReferenceIndex(1);
						setAlignmentStart(200);
						setReadNegativeStrandFlag(true);
						setCigarString("2S2M");
						setReadBases(B("GGGG"));
					}}).getBreakendSummary());
	}
	@Test
	public void incorporateRealignment_should_map_unanchored_breakend() {
		ProcessingContext pc = getContext();
		AssemblyEvidenceSource aes = AES();		
		assertEquals(new BreakpointSummary(0, FWD, 100, 200, 1, BWD, 400, 500),
				AssemblyFactory.incorporateRealignment(pc,
					AssemblyFactory.createUnanchored(pc, aes, Sets.<DirectedEvidence>newHashSet(new MockDirectedEvidence(0, FWD, 100, 200)), B("CCCCC"), B("CCCCC"), 0, 0),
					new SAMRecord(pc.getBasicSamHeader()) {{
						setMappingQuality(40);
						setReferenceIndex(1);
						setAlignmentStart(500);
						setReadNegativeStrandFlag(false);
						setReadBases(B("CCCCC"));
						setCigarString("5M");
					}}).getBreakendSummary());
		
		assertEquals(new BreakpointSummary(0, FWD, 100, 200, 1, FWD, 400, 500),
				AssemblyFactory.incorporateRealignment(pc,
					AssemblyFactory.createUnanchored(pc, aes, Sets.<DirectedEvidence>newHashSet(new MockDirectedEvidence(0, FWD, 100, 200)), B("CCCCC"), B("CCCCC"), 0, 0),
					new SAMRecord(pc.getBasicSamHeader()) {{
						setMappingQuality(40);
						setReferenceIndex(1);
						setAlignmentStart(396);
						setReadNegativeStrandFlag(true);
						setReadBases(B("GGGGG"));
						setCigarString("5M");
					}}).getBreakendSummary());
		
		assertEquals(new BreakpointSummary(0, BWD, 100, 200, 1, FWD, 400, 500),
				AssemblyFactory.incorporateRealignment(pc,
					AssemblyFactory.createUnanchored(pc, aes, Sets.<DirectedEvidence>newHashSet(new MockDirectedEvidence(0, BWD, 100, 200)), B("CCCCC"), B("CCCCC"), 0, 0),
					new SAMRecord(pc.getBasicSamHeader()) {{
						setMappingQuality(40);
						setReferenceIndex(1);
						setAlignmentStart(396);
						setReadNegativeStrandFlag(false);
						setReadBases(B("CCCCC"));
						setCigarString("5M");
					}}).getBreakendSummary());
		
		assertEquals(new BreakpointSummary(0, BWD, 100, 200, 1, BWD, 400, 500),
				AssemblyFactory.incorporateRealignment(pc,
					AssemblyFactory.createUnanchored(pc, aes, Sets.<DirectedEvidence>newHashSet(new MockDirectedEvidence(0, BWD, 100, 200)), B("CCCCC"), B("CCCCC"), 0, 0),
					new SAMRecord(pc.getBasicSamHeader()) {{
						setMappingQuality(40);
						setReferenceIndex(1);
						setAlignmentStart(500);
						setReadNegativeStrandFlag(true);
						setReadBases(B("GGGGG"));
						setCigarString("5M");
					}}).getBreakendSummary());
	}
	@Test
	public void incorporateRealignment_should_not_convert_poorly_mapped_to_breakpoint() {
		SAMRecord ra = Read(1, 100, "1S1M1S");
		ra.setReadBases(B("CGT"));
		ra.setMappingQuality(1);
		ra.setBaseQualities(new byte[] { 0,1,2});
		SAMRecordAssemblyEvidence be = AssemblyFactory.incorporateRealignment(getContext(), big(), ra);
		assertFalse(be instanceof DirectedBreakpoint);
	}
	@Test
	public void getAnchorSequenceString_should_return_entire_assembly_anchor() {
		assertEquals("AAAAT", S(big().getAssemblyAnchorSequence()));
	}
	@Test
	public void breakpoint_should_retain_base_quals() {
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
		SAMRecordAssemblyEvidence e = AssemblyFactory.createAnchored(pc, AES(), BWD, support, 0, 10, 5, B("GGTAAAAC"), new byte[] { 7,6,5,4,3,2,1,0}, 21, 23);
		assertArrayEquals(new byte[] { 7,6,5 }, e.getBreakendQuality());
		SAMRecord ra = Read(1, 102, "1S1M1S");
		ra.setReadBases(B("GGT"));
		ra.setMappingQuality(7);
		ra.setBaseQualities(new byte[] { 0,1,2});
		SAMRecordAssemblyEvidence e2 = AssemblyFactory.incorporateRealignment(getContext(), e, ra);
		assertArrayEquals(e.getBreakendQuality(), e2.getBreakendQuality());
	}
	@Test
	public void getBreakendQual_should_be_sum_of_evidence_quals() {
		double target = 0;
		for (DirectedEvidence e : bigr().getEvidence()) {
			target += e.getBreakendQual();
		}
		assertEquals(target, bigr().getBreakendQual(), 0.001);
	}
	@Test
	public void should_keep_evidence_upon_construction() {
		assertEquals(9, bigr().getEvidence().size());
	}
	@Test
	public void getLocalMapq_should_be_max_evidence_mapq() {
		assertEquals(3, AssemblyFactory.createUnanchored(
				getContext(), AES(), Sets.<DirectedEvidence>newHashSet(
						SCE(FWD, withMapq(1, Read(0, 1, "1M1S"))),
						SCE(FWD, withMapq(2, Read(0, 1, "1M1S"))),
						NRRP(withMapq(3, DP(0, 1, "1M", true, 0, 1, "1M", true))[0], withMapq(4, DP(0, 1, "1M", true, 0, 1, "1M", true))[1])),
						B("T"), B("T"),1, 2).getLocalMapq());
	}
	@Test
	public void id_should_contain_position_direction() {
		Set<DirectedEvidence> evidence = Sets.<DirectedEvidence>newHashSet(
				NRRP(DP(0, 123, "2M", true, 0, 345, "2M", false)));
		AssemblyEvidence e = AssemblyFactory.createUnanchored(
				getContext(), AES(), evidence,
				B("AAAAA"), B("AAAAA"), 2, 0);
		assertTrue(e.getEvidenceID().contains("polyA"));
		assertTrue(e.getEvidenceID().contains(Integer.toString(e.getBreakendSummary().start)));
		assertTrue(e.getEvidenceID().contains(Integer.toString(e.getBreakendSummary().end)));
		assertTrue(e.getEvidenceID().contains("f"));
	}
	@Test
	public void id_be_assembly_unique() {
		AssemblyEvidence e1a = AssemblyFactory.createAnchored(getContext(), AES(), FWD, null, 0, 1, 1, B("AAAAA"), B("AAAAA"), 0, 0);
		AssemblyEvidence e1b = AssemblyFactory.createAnchored(getContext(), AES(), FWD, null, 0, 1, 1, B("AAAAA"), B("AAAAA"), 0, 0);
		assertEquals(e1a.getEvidenceID(), e1b.getEvidenceID());
		for (AssemblyEvidence e : new AssemblyEvidence[] {
				AssemblyFactory.createAnchored(getContext(), AES(), BWD, null, 0, 1, 1, B("AAAAA"), B("AAAAA"), 0, 0),
				AssemblyFactory.createAnchored(getContext(), AES(), FWD, null, 1, 1, 1, B("AAAAA"), B("AAAAA"), 0, 0),
				AssemblyFactory.createAnchored(getContext(), AES(), FWD, null, 0, 2, 1, B("AAAAA"), B("AAAAA"), 0, 0),
				AssemblyFactory.createAnchored(getContext(), AES(), FWD, null, 0, 1, 3, B("AAAAA"), B("AAAAA"), 0, 0),
				AssemblyFactory.createAnchored(getContext(), AES(), FWD, null, 0, 1, 1, B("AAAAT"), B("AAAAA"), 0, 0),
		}) {
			assertNotEquals(e1a.getEvidenceID(), e.getEvidenceID());
		}
	}
	@Test
	public void should_include_untemplated_sequence_for_imprecise_breakpoint() {
		ProcessingContext pc = getContext();
		AssemblyEvidenceSource aes = AES();
		assertEquals("GT", 
				((RealignedSAMRecordAssemblyEvidence)AssemblyFactory.incorporateRealignment(pc,
					AssemblyFactory.createUnanchored(pc, aes, Sets.<DirectedEvidence>newHashSet(new MockDirectedEvidence(0, FWD, 100, 200)), B("GTNAC"), B("CCCCC"), 0, 0),
					new SAMRecord(pc.getBasicSamHeader()) {{
						setMappingQuality(40);
						setReferenceIndex(1);
						setAlignmentStart(500);
						setReadNegativeStrandFlag(false);
						setReadBases(B("GTNAC"));
						setCigarString("2S3M");
					}})).getUntemplatedSequence());
	}
	@Test
	public void should_set_breakend_exact() {
		assertTrue(AssemblyFactory.createAnchored(getContext(), AES(), FWD, null, 0, 1, 1, B("AAAAA"), B("AAAAA"), 0, 0).isBreakendExact());
		assertFalse(AssemblyFactory.createUnanchored(getContext(), AES(), Sets.<DirectedEvidence>newHashSet(NRRP(DP(0, 10, "2M", true, 0, 10, "2M", false))), B("AAAAA"), B("AAAAA"), 2, 0)
				.isBreakendExact());
	}
}
