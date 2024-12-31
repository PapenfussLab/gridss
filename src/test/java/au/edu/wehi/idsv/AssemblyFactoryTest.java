package au.edu.wehi.idsv;

import au.edu.wehi.idsv.sam.SAMRecordUtil;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Lists;
import htsjdk.samtools.SAMRecord;
import org.junit.Test;

import java.util.List;
import java.util.stream.DoubleStream;

import static org.junit.Assert.*;


public class AssemblyFactoryTest extends TestHelper {
	@Test
	public void should_set_breakend_anchored_bwd() {
		SingleReadEvidence e = asAssemblyEvidence(AssemblyFactory.createAnchoredBreakend(
				getContext(), AES(), new SequentialIdGenerator("asm"), BWD, null,null,
				0, 1, 2, B("GTACCC"), new byte[] { 1, 2, 3, 4, 4, 8 }));
		assertEquals(new BreakendSummary(0, BWD, 1), e.getBreakendSummary());
	}
	@Test
	public void should_set_breakend_anchored_fwd() {
		SingleReadEvidence e = asAssemblyEvidence(AssemblyFactory.createAnchoredBreakend(
				getContext(), AES(), new SequentialIdGenerator("asm"), FWD, null,null,
				0, 10, 2, B("GTACCC"), new byte[] { 1, 2, 3, 4, 4, 8 }));
		assertEquals(new BreakendSummary(0, FWD, 10), e.getBreakendSummary());
	}
	@Test
	public void should_set_assembly_properties_bwd() {
		SingleReadEvidence e = asAssemblyEvidence(AssemblyFactory.createAnchoredBreakend(
				getContext(), AES(), new SequentialIdGenerator("asm"), BWD, null,null,
				0, 1, 2, B("GTACCC"), new byte[] { 1, 2, 3, 4, 4, 8}));
		assertArrayEquals(new byte[] { 1, 2, 3, 4, }, e.getBreakendQuality());
		assertEquals("GTAC", S(e.getBreakendSequence()));
		assertEquals("CC", S(e.getAnchorSequence()));
		assertEquals("GTACCC", S(e.getSAMRecord().getReadBases()));
	}
	@Test
	public void should_set_assembly_properties_fwd() {
		SingleReadEvidence e = asAssemblyEvidence(AssemblyFactory.createAnchoredBreakend(
				getContext(), AES(), new SequentialIdGenerator("asm"), FWD, null,null,
				0, 2, 2, B("GTACCC"), new byte[] { 1, 3, 3, 4, 4, 8}));
		assertArrayEquals(new byte[] { 3, 4, 4, 8, }, e.getBreakendQuality());
		assertEquals("ACCC", S(e.getBreakendSequence()));
		assertEquals("GT", S(e.getAnchorSequence()));
		assertEquals("GTACCC", S(e.getSAMRecord().getReadBases()));
	}
	@Test
	public void truncate_should_removed_anchor_bases() {
		SingleReadEvidence e = asAssemblyEvidence(AssemblyFactory.createAnchoredBreakend(
				getContext(), AES(), new SequentialIdGenerator("asm"), FWD, null,null,
				0, 1, 2, B("GTACCC"), new byte[] { 1, 3, 3, 4, 4, 8}));
		assertEquals(1, e.getSAMRecord().getAlignmentStart());
		assertEquals("1M4S", e.getSAMRecord().getCigarString());
		assertArrayEquals(new byte[] { 3, 4, 4, 8, }, e.getBreakendQuality());
		assertEquals("ACCC", S(e.getBreakendSequence()));
		assertEquals("T", S(e.getAnchorSequence()));
		assertEquals("TACCC", S(e.getSAMRecord().getReadBases()));
	}
	@Test
	public void should_not_truncate_if_truncation_removes_cigar_element() {
		SingleReadEvidence e = asAssemblyEvidence(AssemblyFactory.createAnchoredBreakend(
				getContext(), AES(), new SequentialIdGenerator("asm"), FWD, null,null,
				0, 100000, 2, B("GTACCC"), new byte[] { 1, 3, 3, 4, 4, 8}));
		assertEquals(99999, e.getSAMRecord().getAlignmentStart());
		assertEquals("2M4S", e.getSAMRecord().getCigarString());
	}
	@Test
	public void should_set_breakend_unanchored() {
		for (BreakendSummary bs : ImmutableList.of(
				new BreakendSummary(0, FWD, 2, 1, 2),
				new BreakendSummary(0, BWD, 1, 1, 2)
				)) {
			SingleReadEvidence e = asAssemblyEvidence(AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), bs, null, null, B("AAAAA"), B("AAAAA") ));
			assertEquals(bs, e.getBreakendSummary());
		}
	}
	@Test
	public void createUnanchoredBreakend_should_round_trip_unanchored_breakend_start_end_but_not_nominal_position() {
		assertEquals(new BreakendSummary(0, FWD, 2, 1, 3), asAssemblyEvidence(AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), new BreakendSummary(0, FWD, 1, 1, 3), null, null, B("GTAC"), new byte[] {1,2,3,4})).getBreakendSummary());
		assertEquals(new BreakendSummary(0, BWD, 15, 10, 20), asAssemblyEvidence(AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), new BreakendSummary(0, BWD, 11, 10, 20), null, null, B("GTAC"), new byte[] {1,2,3,4})).getBreakendSummary());
	}

	@Test
	public void should_set_assembly_attribute_ASSEMBLY_MAPQ_REMOTE_MAX() {
		assertEquals(17, bigr().getRemoteMapq());
	}
	@Test
	public void should_set_assembly_attribute_ASSEMBLY_LENGTH_LOCAL_MAX() {
		assertEquals(5, big().getAnchorSequence().length);
		assertEquals(5, bigr().getAnchorSequence().length);
	}
	@Test
	public void should_set_assembly_attribute_ASSEMBLY_LENGTH_REMOTE_MAX() {
		assertEquals(3, big().getBreakendSequence().length);
		assertEquals(3, bigr().getBreakendSequence().length);
	}
	@Test
	public void should_set_assembly_attribute_ASSEMBLY_READPAIR_COUNT() {
		assertEquals(2, new AssemblyAttributes(bigr()).getSupportingReadCount(bigr().getBreakendAssemblyContigOffset(), ImmutableSet.of(0), ImmutableSet.of(AssemblyEvidenceSupport.SupportType.ReadPair), (AssemblyEvidenceSource)bigr().source, null));
		assertEquals(4, new AssemblyAttributes(bigr()).getSupportingReadCount(bigr().getBreakendAssemblyContigOffset(), ImmutableSet.of(1), ImmutableSet.of(AssemblyEvidenceSupport.SupportType.ReadPair), (AssemblyEvidenceSource)bigr().source, null));
	}
	@Test
	public void should_set_assembly_attribute_ASSEMBLY_READPAIR_LENGTH_MAX() {
		assertEquals(10, new AssemblyAttributes(bigr()).getAssemblyMaxReadLength());
	}
	@Test
	public void should_set_assembly_attribute_ASSEMBLY_SOFTCLIP_COUNT() {
		assertEquals(1, new AssemblyAttributes(bigr()).getSupportingReadCount(bigr().getBreakendAssemblyContigOffset(), ImmutableSet.of(0), ImmutableSet.of(AssemblyEvidenceSupport.SupportType.Read), (AssemblyEvidenceSource)bigr().source, null));
		assertEquals(2, new AssemblyAttributes(bigr()).getSupportingReadCount(bigr().getBreakendAssemblyContigOffset(), ImmutableSet.of(1), ImmutableSet.of(AssemblyEvidenceSupport.SupportType.Read), (AssemblyEvidenceSource)bigr().source, null));
	}
//	@Test
//	public void should_set_assembly_attribute_ASSEMBLY_REMOTE_COUNT() {
//		ProcessingContext pc = getContext();
//		MockSAMEvidenceSource tes = new MockSAMEvidenceSource(pc);
//		tes.category = 1;
//		Set<DirectedEvidence> support = Sets.newHashSet();
//		support.add(SR(Read(1, 10, "4S1M"), Read(0, 10, "4M1S")));
//		SAMRecord ass = AssemblyFactory.createAnchoredBreakend(pc, AES(), BWD, support, 0, 10, 5, B("CGTAAAAT"), new byte[] { 0,1,2,3,4,5,6,7});
//		assertEquals(0, new AssemblyAttributes(ass).getAssemblySupportCountSoftClipRemote(0));
//		assertEquals(1, new AssemblyAttributes(ass).getAssemblySupportCountSoftClipRemote(1));
//	}
	@Test
	public void should_set_assembly_getLocalMapq() {
		assertEquals(10, bigr().getLocalMapq());
	}
	@Test
	public void should_set_assembly_attribute_ASSEMBLY_BREAKEND_QUALS() {
		assertArrayEquals(new byte[] { 0,1,2}, big().getBreakendQuality());
	}
	@Test
	public void should_set_evidence_source() {
		AssemblyEvidenceSource es = AES();
		SingleReadEvidence e = asAssemblyEvidence(es, AssemblyFactory.createAnchoredBreakend(
				es.getContext(), es, new SequentialIdGenerator("asm"), FWD, null, null,
				0, 1, 2, B("GTACCC"), new byte[] { 1, 3, 3, 4, 4, 8}));
		assertEquals(es, e.getEvidenceSource());
	}
	public SoftClipEvidence big() {
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
		SAMRecord ass = AssemblyFactory.createAnchoredBreakend(pc, AES(), new SequentialIdGenerator("asm"), BWD, support, null, 0, 10, 5, B("CGTAAAAT"), new byte[] { 0,1,2,3,4,5,6,7});
		return SoftClipEvidence.create(AES(), BWD, ass);
	}
	public SplitReadEvidence bigr() {
		SAMRecord r = big().getSAMRecord();
		SAMRecord ra = Read(1, 100, "1S1M1S5S");
		ra.setReadBases(B("CGTAAAAT"));
		ra.setMappingQuality(17);
		ra.setBaseQualities(new byte[] { 0,1,2,3,4,5,6,7 });
		ra.setReadName(SplitReadHelper.getSplitReadRealignments(r, false, AES().getContext().getEvidenceIDGenerator(), (byte)0).get(0).getReadName());
		SplitReadHelper.convertToSplitRead(r, ImmutableList.of(ra), null, false);
		return SplitReadEvidence.create(AES(), r).get(0);
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
		
		assertEquals(new BreakpointSummary(0, FWD, 100, 1, BWD, 200),
				incorporateRealignment(aes,
					AssemblyFactory.createAnchoredBreakend(pc, aes, new SequentialIdGenerator("asm"), FWD, null, null, 0, 100, 1, B("NNN"), B("NNN")),
					ImmutableList.of(new SAMRecord(pc.getBasicSamHeader()) {/**
						 * 
						 */
						private static final long serialVersionUID = 1L;

					{
						setMappingQuality(40);
						setReferenceIndex(1);
						setAlignmentStart(200);
						setReadNegativeStrandFlag(false);
						setCigarString("2M");
						setReadBases(B("NN"));
					}})).getBreakendSummary());
		
		assertEquals(new BreakpointSummary(0, FWD, 100, 1, FWD, 200),
				incorporateRealignment(aes,
					AssemblyFactory.createAnchoredBreakend(pc, aes, new SequentialIdGenerator("asm"), FWD, null, null, 0, 100, 1, B("NNN"), B("NNN")),
					ImmutableList.of(new SAMRecord(pc.getBasicSamHeader()) {/**
						 * 
						 */
						private static final long serialVersionUID = 1L;

					{
						setMappingQuality(40);
						setReferenceIndex(1);
						setAlignmentStart(199);
						setReadNegativeStrandFlag(true);
						setCigarString("2M");
						setReadBases(B("NN"));
					}})).getBreakendSummary());
		
		assertEquals(new BreakpointSummary(0, BWD, 100, 1, FWD, 200),
				incorporateRealignment(aes,
					AssemblyFactory.createAnchoredBreakend(pc, aes, new SequentialIdGenerator("asm"), BWD, null, null, 0, 100, 1, B("NNN"), B("NNN")),
					ImmutableList.of(new SAMRecord(pc.getBasicSamHeader()) {/**
						 * 
						 */
						private static final long serialVersionUID = 1L;

					{
						setMappingQuality(40);
						setReferenceIndex(1);
						setAlignmentStart(199);
						setReadNegativeStrandFlag(false);
						setCigarString("2M");
						setReadBases(B("NN"));
					}})).getBreakendSummary());
		
		assertEquals(new BreakpointSummary(0, BWD, 100, 1, BWD, 200),
				incorporateRealignment(aes,
					AssemblyFactory.createAnchoredBreakend(pc, aes, new SequentialIdGenerator("asm"), BWD, null, null, 0, 100, 1, B("NNN"), B("NNN")),
					ImmutableList.of(new SAMRecord(pc.getBasicSamHeader()) {/**
						 * 
						 */
						private static final long serialVersionUID = 1L;

					{
						setMappingQuality(40);
						setReferenceIndex(1);
						setAlignmentStart(200);
						setReadNegativeStrandFlag(true);
						setCigarString("2M");
						setReadBases(B("NN"));
					}})).getBreakendSummary());
	}
	@Test
	public void incorporateRealignment_should_map_anchored_breakend_with_soft_clip() {
		ProcessingContext pc = getContext();
		AssemblyEvidenceSource aes = AES();
		
		assertEquals(new BreakpointSummary(0, FWD, 100, 1, BWD, 200),
				incorporateRealignment(aes,
					AssemblyFactory.createAnchoredBreakend(pc, aes, new SequentialIdGenerator("asm"), FWD, null, null, 0, 100, 1, B("CCCCC"), B("CCCCC")),
					ImmutableList.of(new SAMRecord(pc.getBasicSamHeader()) {/**
						 * 
						 */
						private static final long serialVersionUID = 1L;

					{
						setMappingQuality(40);
						setReferenceIndex(1);
						setAlignmentStart(200);
						setReadNegativeStrandFlag(false);
						setReadBases(B("CCCC"));
						setCigarString("2S2M");
					}})).getBreakendSummary());
		
		assertEquals(new BreakpointSummary(0, FWD, 100, 1, FWD, 200),
				incorporateRealignment(aes,
					AssemblyFactory.createAnchoredBreakend(pc, aes, new SequentialIdGenerator("asm"), FWD, null, null, 0, 100, 1, B("CCCCC"), B("CCCCC")),
					ImmutableList.of(new SAMRecord(pc.getBasicSamHeader()) {/**
						 * 
						 */
						private static final long serialVersionUID = 1L;

					{
						setMappingQuality(40);
						setReferenceIndex(1);
						setAlignmentStart(199);
						setReadNegativeStrandFlag(true);
						setReadBases(B("GGGG"));
						setCigarString("2M2S");
					}})).getBreakendSummary());
		
		assertEquals(new BreakpointSummary(0, BWD, 100, 1, FWD, 200),
				incorporateRealignment(aes,
					AssemblyFactory.createAnchoredBreakend(pc, aes, new SequentialIdGenerator("asm"), BWD, null, null, 0, 100, 1, B("CCCCC"), B("CCCCC")),
					ImmutableList.of(new SAMRecord(pc.getBasicSamHeader()) {/**
						 * 
						 */
						private static final long serialVersionUID = 1L;

					{
						setMappingQuality(40);
						setReferenceIndex(1);
						setAlignmentStart(199);
						setReadNegativeStrandFlag(false);
						setReadBases(B("CCCC"));
						setReadBases(B("GGGG"));
						setCigarString("2M2S");
					}})).getBreakendSummary());
		
		assertEquals(new BreakpointSummary(0, BWD, 100, 1, BWD, 200),
				incorporateRealignment(aes,
					AssemblyFactory.createAnchoredBreakend(pc, aes, new SequentialIdGenerator("asm"), BWD, null, null, 0, 100, 1, B("CCCCC"), B("CCCCC")),
					ImmutableList.of(new SAMRecord(pc.getBasicSamHeader()) {/**
						 * 
						 */
						private static final long serialVersionUID = 1L;

					{
						setMappingQuality(40);
						setReferenceIndex(1);
						setAlignmentStart(200);
						setReadNegativeStrandFlag(true);
						setCigarString("2S2M");
						setReadBases(B("GGGG"));
					}})).getBreakendSummary());
	}
	@Test
	public void incorporateRealignment_should_map_unanchored_breakend() {
		ProcessingContext pc = getContext();
		AssemblyEvidenceSource aes = AES();		
		assertEquals(new BreakpointSummary(0, FWD, 150, 100, 200, 1, BWD, 450, 400, 500),
				incorporateRealignment(aes,
					AssemblyFactory.createUnanchoredBreakend(pc, aes, new SequentialIdGenerator("asm"), new BreakendSummary(0, FWD, 150, 100, 200), null, null, B("CCCCC"), B("CCCCC")),
					ImmutableList.of(new SAMRecord(pc.getBasicSamHeader()) {/***/
						private static final long serialVersionUID = 1L;
					{
						setMappingQuality(40);
						setReferenceIndex(1);
						setAlignmentStart(500);
						setReadNegativeStrandFlag(false);
						setReadBases(B("CCCCC"));
						setCigarString("5M");
					}})).getBreakendSummary());
		
		assertEquals(new BreakpointSummary(0, FWD, 150, 100, 200, 1, FWD, 450, 400, 500),
				incorporateRealignment(aes,
					AssemblyFactory.createUnanchoredBreakend(pc, aes, new SequentialIdGenerator("asm"), new BreakendSummary(0, FWD, 150, 100, 200), null, null, B("CCCCC"), B("CCCCC")),
					ImmutableList.of(new SAMRecord(pc.getBasicSamHeader()) {/**
						 * 
						 */
						private static final long serialVersionUID = 1L;

					{
						setMappingQuality(40);
						setReferenceIndex(1);
						setAlignmentStart(396);
						setReadNegativeStrandFlag(true);
						setReadBases(B("GGGGG"));
						setCigarString("5M");
					}})).getBreakendSummary());
		
		assertEquals(new BreakpointSummary(0, BWD, 150, 100, 200, 1, FWD, 450, 400, 500),
				incorporateRealignment(aes,
					AssemblyFactory.createUnanchoredBreakend(pc, aes, new SequentialIdGenerator("asm"), new BreakendSummary(0, BWD, 150, 100, 200), null, null, B("CCCCC"), B("CCCCC")),
					ImmutableList.of(new SAMRecord(pc.getBasicSamHeader()) {/**
						 * 
						 */
						private static final long serialVersionUID = 1L;

					{
						setMappingQuality(40);
						setReferenceIndex(1);
						setAlignmentStart(396);
						setReadNegativeStrandFlag(false);
						setReadBases(B("CCCCC"));
						setCigarString("5M");
					}})).getBreakendSummary());
		
		assertEquals(new BreakpointSummary(0, BWD, 150, 100, 200, 1, BWD, 450, 400, 500),
				incorporateRealignment(aes,
					AssemblyFactory.createUnanchoredBreakend(pc, aes, new SequentialIdGenerator("asm"), new BreakendSummary(0, BWD, 150, 100, 200), null, null, B("CCCCC"), B("CCCCC")),
					ImmutableList.of(new SAMRecord(pc.getBasicSamHeader()) {/**
						 * 
						 */
						private static final long serialVersionUID = 1L;

					{
						setMappingQuality(40);
						setReferenceIndex(1);
						setAlignmentStart(500);
						setReadNegativeStrandFlag(true);
						setReadBases(B("GGGGG"));
						setCigarString("5M");
					}})).getBreakendSummary());
	}
	@Test
	public void getAnchorSequenceString_should_return_entire_assembly_anchor() {
		assertEquals("AAAAT", S(big().getAnchorSequence()));
	}
	@Test
	public void breakpoint_should_retain_base_quals() {
		ProcessingContext pc = getContext();
		SingleReadEvidence e = asAssemblyEvidence(AssemblyFactory.createAnchoredBreakend(pc, AES(), new SequentialIdGenerator("asm"), BWD, null, null, 0, 10, 5, B("GGTAAAAC"), new byte[] { 7,6,5,4,3,2,1,0}));
		assertArrayEquals(new byte[] { 7,6,5 }, e.getBreakendQuality());
		SAMRecord ra = Read(1, 102, "1S1M1S");
		ra.setReadBases(B("GGT"));
		ra.setMappingQuality(7);
		ra.setBaseQualities(new byte[] { 0,1,2});
		SingleReadEvidence e2 = incorporateRealignment(AES(), e, ImmutableList.of(ra));
		assertArrayEquals(e.getBreakendQuality(), e2.getBreakendQuality());
	}
	@Test
	public void getQual_should_be_sum_of_evidence() {
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
		SAMRecord ass = AssemblyFactory.createAnchoredBreakend(pc, AES(), new SequentialIdGenerator("asm"), BWD, support, fullSupport(support), 0, 10, 5, B("CGTAAAAT"), new byte[] { 0,1,2,3,4,5,6,7});
		SoftClipEvidence asse = SoftClipEvidence.create(AES(), BWD, ass);

		double[] supportQual = support.stream().mapToDouble(e -> e.getBreakendQual()).toArray();
		assertEquals(DoubleStream.of(supportQual).sum(), asse.getBreakendQual(), 0.001);
	}
	@Test
	public void getLocalMapq_should_be_max_evidence_mapq() {
		List<DirectedEvidence> support = Lists.<DirectedEvidence>newArrayList(
				SCE(FWD, withMapq(10, Read(0, 1, "1M2S"))),
				SCE(FWD, withMapq(20, Read(0, 1, "1M1S"))),
				NRRP(withMapq(30, DP(0, 1, "1M", true, 1, 1, "1M", true))[0], withMapq(4, DP(0, 1, "1M", true, 1, 1, "1M", true))[1]));
		SingleReadEvidence ass = asAssemblyEvidence(AssemblyFactory.createUnanchoredBreakend(
				getContext(), AES(), new SequentialIdGenerator("asm"), new BreakendSummary(0, FWD, 1), support, null, B("T"), B("T")));
		assertEquals(30, ass.getLocalMapq());
	}
	@Test
	public void id_be_assembly_unique() {
		ProcessingContext context = getContext();
		AssemblyEvidenceSource aes = AES();
		SequentialIdGenerator gen = new SequentialIdGenerator("asm");
		SingleReadEvidence e1a = asAssemblyEvidence(AssemblyFactory.createAnchoredBreakend(context, aes, gen, FWD, null, null, 0, 1, 1, B("AAAAA"), B("AAAAA")));
		SingleReadEvidence e1b = asAssemblyEvidence(AssemblyFactory.createAnchoredBreakend(context, aes, gen, FWD, null, null, 0, 1, 1, B("AAAAA"), B("AAAAA")));
		assertNotEquals(e1a.getEvidenceID(), e1b.getEvidenceID());
		for (SingleReadEvidence e : new SingleReadEvidence[] {
				asAssemblyEvidence(AssemblyFactory.createAnchoredBreakend(context, aes, gen, BWD, null, null, 0, 1, 1, B("AAAAA"), B("AAAAA"))),
				asAssemblyEvidence(AssemblyFactory.createAnchoredBreakend(context, aes, gen, FWD, null, null, 1, 1, 1, B("AAAAA"), B("AAAAA"))),
				asAssemblyEvidence(AssemblyFactory.createAnchoredBreakend(context, aes, gen, FWD, null, null, 0, 2, 1, B("AAAAA"), B("AAAAA"))),
				asAssemblyEvidence(AssemblyFactory.createAnchoredBreakend(context, aes, gen, FWD, null, null, 0, 1, 3, B("AAAAA"), B("AAAAA"))),
				asAssemblyEvidence(AssemblyFactory.createAnchoredBreakend(context, aes, gen, FWD, null, null, 0, 1, 1, B("AAAAT"), B("AAAAA"))),
		}) {
			assertNotEquals(e1a.getEvidenceID(), e.getEvidenceID());
		}
	}
	@Test
	public void should_include_untemplated_sequence_for_imprecise_breakpoint() {
		ProcessingContext pc = getContext();
		AssemblyEvidenceSource aes = AES();
		assertEquals("GT", 
				incorporateRealignment(aes,
					AssemblyFactory.createUnanchoredBreakend(pc, aes, new SequentialIdGenerator("asm"), new BreakendSummary(0, FWD, 150, 100, 200), null, null, B("GTNAC"), B("CCCCC")),
					ImmutableList.of(new SAMRecord(pc.getBasicSamHeader()) {/**
						 * 
						 */
						private static final long serialVersionUID = 1L;

					{
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
		assertTrue(asAssemblyEvidence(AssemblyFactory.createAnchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), FWD, null, null, 0, 1, 1, B("AAAAA"), B("AAAAA"))).isBreakendExact());
		assertFalse(asAssemblyEvidence(AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), new BreakendSummary(0, FWD, 150, 100, 200), null, null, B("AAAAA"), B("AAAAA"))).isBreakendExact());
	}
	@Test
	public void should_assemble_breakpoint() {
		SAMRecord r = AssemblyFactory.createAnchoredBreakpoint(getContext(), AES(), new SequentialIdGenerator("asm"),
				null,null,
				0, 1, 1, 0, 2, 1,
				B("NAAAN"), B("AAAAA"));
		IndelEvidence e = IndelEvidence.create(AES(), 0, r).get(0);
		assertTrue(e.getBreakendSummary() instanceof BreakpointSummary);
		BreakpointSummary bp = (BreakpointSummary)e.getBreakendSummary();
		assertEquals(0, bp.referenceIndex);
		assertEquals(1, bp.start);
		assertEquals(1, bp.end);
		assertEquals(FWD, bp.direction);
		assertEquals(0, bp.referenceIndex2);
		assertEquals(2, bp.start2);
		assertEquals(2, bp.end2);
		assertEquals(BWD, bp.direction2);
		assertEquals("AAA", e.getUntemplatedSequence());
	}
	@Test
	public void breakpoint_assembly_should_convert_tandem_duplication_assembly_soft_clip_to_ensure_cigar_is_valid() {
		SAMRecord r = AssemblyFactory.createAnchoredBreakpoint(getContext(), AES(), new SequentialIdGenerator("asm"), null,null,
				0, 100, 1,
				0, 1, 2, B("AAA"), B("AAA"));
		List<SingleReadEvidence> e = SingleReadEvidence.createEvidence(AES(), 0, r);
		assertEquals(1, e.size());
		assertTrue(e.get(0) instanceof SoftClipEvidence);
	}
	@Test
	public void createAnchoredBreakpoint_should_not_crash_on_ref_allele_breakpoint() {
		// 1M1M
		SAMRecord r = AssemblyFactory.createAnchoredBreakpoint(getContext(), AES(), new SequentialIdGenerator("asm"), null, null, 0, 10, 1, 0, 11, 1, B("AA"), B("AA"));
		assertTrue(SAMRecordUtil.isReferenceAlignment(r));
	}
	@Test
	public void imprecise_should_include_ci_interval_on_both_sides() {
		SingleReadEvidence e = asAssemblyEvidence(AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), new BreakendSummary(0, BWD, 1, 1, 10), null, null, B("AA"), B("AA")));
		SplitReadEvidence re = (SplitReadEvidence)incorporateRealignment(AES(), e.getSAMRecord(), ImmutableList.of(Read(1, 100, "2M")));
		assertEquals(new BreakpointSummary(0, BWD, 5, 1, 10, 1, FWD, 101, 101, 110), re.getBreakendSummary());
	}
	@Test
	public void realignment_to_negative_strand_should_change_anchor_position() {
		SingleReadEvidence e = asAssemblyEvidence(AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), new BreakendSummary(0, BWD, 1, 1, 10), null, null, B("AAAAAAAAAA"), B("0000000000")));
		SplitReadEvidence re = (SplitReadEvidence)incorporateRealignment(AES(), e.getSAMRecord(), ImmutableList.of(onNegative(Read(1, 200, "10M"))[0]));
		assertEquals(new BreakpointSummary(0, BWD, 5, 1, 10, 1, BWD, 195, 191, 200), re.getBreakendSummary());
	}
	@Test
	public void realignment_of_breakpoint_outside_of_contig_bounds_should_have_inexact_breakpoint() {
		SAMRecord r = AssemblyFactory.createUnanchoredBreakend(
				getContext(), AES(), new SequentialIdGenerator("asm"), new BreakendSummary(0, BWD, -20, -30, -10),
				null, null, B("GTA"), new byte[] { 1, 3, 3});
		// 20
		//  012345678
		//  MMM
		//    |----  .. ---| 
		SplitReadEvidence re = (SplitReadEvidence)incorporateRealignment(AES(), r, ImmutableList.of(Read(1, 200, "3M")));
		assertEquals(new BreakendSummary(0, BWD, 1), re.getBreakendSummary().localBreakend());
		assertEquals(new BreakendSummary(1, FWD, 212, 202, 222), re.getBreakendSummary().remoteBreakend());
	}
	@Test
	public void basequal_range_should_be_limited() {
		SAMRecord r = AssemblyFactory.createAnchoredBreakend(
				getContext(), AES(), new SequentialIdGenerator("asm"), BWD, null,null,
				0, 1, 2, B("GTACC"), new byte[] { 92,93,(byte)254,-1,94 });
		assertArrayEquals(new byte[] { 92, 93,93,93,93}, r.getBaseQualities());
	}
}

