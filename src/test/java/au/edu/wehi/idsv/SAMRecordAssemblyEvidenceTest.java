package au.edu.wehi.idsv;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;

import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Range;
import org.apache.commons.configuration.ConfigurationException;
import org.junit.Ignore;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import com.google.common.collect.Lists;

import au.edu.wehi.idsv.configuration.GridssConfiguration;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.sam.SamTags;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.Header;


public class SAMRecordAssemblyEvidenceTest extends TestHelper {
	@Test
	public void should_create_SAMRecord_for_assembly() {
		SAMRecord e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), FWD, null, null,
			1, 2, 1, B("GTAC"), new byte[] {1,2,3,4});
		assertNotNull(e);
		assertEquals("GTAC", S(e.getReadBases()));
		assertArrayEquals( new byte[] {1,2,3,4}, e.getBaseQualities());
	}
	@Test
	public void anchor_positions_should_match_genomic() {
		SAMRecord fwd = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), FWD, null, null, 1, 10, 3, B("GTACCCA"), new byte[] { 1, 2, 3, 4, 4, 8, 8 });
		SAMRecord bwd = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), BWD, null, null, 1, 10, 3, B("GTACCCA"), new byte[] { 1, 2, 3, 4, 4, 8, 8 });
		assertEquals(1, (int)fwd.getReferenceIndex());
		assertEquals(1, (int)bwd.getReferenceIndex());
		assertEquals(8, fwd.getAlignmentStart());
		assertEquals(10, bwd.getAlignmentStart());
		assertEquals("3M4S", fwd.getCigarString());
		assertEquals("4S3M", bwd.getCigarString());
		assertEquals("CCCA", S(asAssemblyEvidence(fwd).getBreakendSequence()));
		assertEquals("GTAC", S(asAssemblyEvidence(bwd).getBreakendSequence()));
		assertEquals("GTA", S(asAssemblyEvidence(fwd).getAnchorSequence()));
		assertEquals("CCA", S(asAssemblyEvidence(bwd).getAnchorSequence()));
	}
	@Test
	public void unanchor_positions_should_match_genomic() {
		SAMRecord fwd = AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), new BreakendSummary(1, FWD, 7, 5, 10), null, null, B("AAA"), B("AAA"));
		SAMRecord bwd = AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), new BreakendSummary(1, BWD, 7, 5, 10), null, null, B("AAA"), B("AAA"));
		assertEquals(1, (int)fwd.getReferenceIndex());
		assertEquals(1, (int)bwd.getReferenceIndex());
		assertEquals(5, fwd.getAlignmentStart());
		assertEquals(5, bwd.getAlignmentStart());
		assertEquals(10, fwd.getAlignmentEnd());
		assertEquals(10, bwd.getAlignmentEnd());
		assertEquals("1X4N1X3S", fwd.getCigarString());
		assertEquals("3S1X4N1X", bwd.getCigarString());
		assertEquals("AAA", S(asAssemblyEvidence(fwd).getBreakendSequence()));
		assertEquals("AAA", S(asAssemblyEvidence(bwd).getBreakendSequence()));
		assertEquals("", S(asAssemblyEvidence(fwd).getAnchorSequence()));
		assertEquals("", S(asAssemblyEvidence(bwd).getAnchorSequence()));
	}
	@Test
	public void unanchor_sequences_should_match_assembly() {
		for (SAMRecord e : new SAMRecord[] {
				AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), new BreakendSummary(1, FWD, 1, 1, 1), null, null, B("AAA"), B("AAA")),
				AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), new BreakendSummary(1, FWD, 1, 1, 2), null, null, B("AAA"), B("AAA")),
				AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), new BreakendSummary(1, FWD, 1, 1, 3), null, null, B("AAA"), B("AAA")),
				AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), new BreakendSummary(1, FWD, 1, 1, 4), null, null, B("AAA"), B("AAA")),
				AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), new BreakendSummary(1, BWD, 1, 1, 1), null, null, B("AAA"), B("AAA")),
				AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), new BreakendSummary(1, BWD, 1, 1, 2), null, null, B("AAA"), B("AAA")),
				AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), new BreakendSummary(1, BWD, 1, 1, 3), null, null, B("AAA"), B("AAA")),
				AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), new BreakendSummary(1, BWD, 1, 1, 4), null, null, B("AAA"), B("AAA"))
		}) {
			assertEquals("AAA", S(asAssemblyEvidence(e).getBreakendSequence()));
			assertEquals("", S(asAssemblyEvidence(e).getAnchorSequence()));
		}
	}
	private void assertMateFields(SAMRecord r, SAMRecord mate) {
		assertEquals(r.getMateNegativeStrandFlag(), mate.getReadNegativeStrandFlag());
		assertEquals(r.getMateUnmappedFlag(), mate.getReadUnmappedFlag());
		assertEquals(r.getMateReferenceIndex(), mate.getReferenceIndex());
		assertEquals(r.getMateAlignmentStart(), mate.getAlignmentStart());
	}
	@Test
	public void should_use_MS_for_anchored_breakend() {
		assertEquals("1M3S", AssemblyFactory.createAnchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), FWD, null, null,
				1, 2, 1, B("GTAC"), new byte[] {1,2,3,4}).getCigarString());
		assertEquals("3S1M", AssemblyFactory.createAnchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), BWD, null, null,
				1, 2, 1, B("GTAC"), new byte[] {1,2,3,4}).getCigarString());
	}
	@Test
	public void should_use_XNXS_for_unanchored_interval_breakend() {
		assertEquals("1X1N1X4S", AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), new BreakendSummary(0, FWD, 1, 1, 3), null, null, B("GTAC"), new byte[] {1,2,3,4}).getCigarString());
		assertEquals("4S1X9N1X", AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), new BreakendSummary(0, BWD, 10, 10, 20), null, null, B("GTAC"), new byte[] {1,2,3,4}).getCigarString());
	}
	@Test
	public void should_use_minimal_cigar_representation() {
		assertEquals("1X3S", AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), new BreakendSummary(1, FWD, 5, 5, 5), null, null, B("AAA"), B("AAA")).getCigarString());
		assertEquals("2X3S", AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), new BreakendSummary(1, FWD, 5, 5, 6), null, null, B("AAA"), B("AAA")).getCigarString());
		assertEquals("1X1N1X3S", AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), new BreakendSummary(1, FWD, 5, 5, 7), null, null, B("AAA"), B("AAA")).getCigarString());
	}
	@Test
	public void should_track_breakend_evidence() {
		DirectedEvidence e1 = SCE(BWD, Read(0, 1, "5S5M"));
		DirectedEvidence e2 = SCE(BWD, Read(0, 1, "6S5M"));
		DirectedEvidence e3 = NRRP(OEA(0, 1, "1M", false));
		DirectedEvidence b1 = SCE(BWD, Read(0, 1, "7S5M"));
		SAMRecord r = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), FWD, Lists.newArrayList(e1, e2, e3), null,
			1, 2, 1, B("GTAC"), new byte[] {1,2,3,4});
		AssemblyAttributes e = new AssemblyAttributes(r);
		assertTrue(e.isPartOfAssembly(e1));
		assertTrue(e.isPartOfAssembly(e2));
		assertTrue(e.isPartOfAssembly(e3));
		assertFalse(e.isPartOfAssembly(b1));
		
		assertEquals(2, (int)e.getSupportingReadCount(1, null, ImmutableSet.of(AssemblyEvidenceSupport.SupportType.Read)));
		assertEquals(1, (int)e.getSupportingReadCount(1, null, ImmutableSet.of(AssemblyEvidenceSupport.SupportType.ReadPair)));
	}
	@Test
	public void should_track_read_name() {
		DirectedEvidence e1 = SCE(BWD, withReadName("r1", Read(0, 1, "5S5M")));
		DirectedEvidence e2 = SCE(BWD, withReadName("r2", Read(0, 1, "6S5M")));
		DirectedEvidence e3 = NRRP(withReadName("r3", OEA(0, 1, "1M", false)));
		MockSAMEvidenceSource cat2 = SES();
		cat2.category = 1;
		DirectedEvidence b1 = SCE(BWD, cat2, withReadName("r4", Read(0, 1, "7S5M")));
		SAMRecord r = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), FWD, Lists.newArrayList(e1, e2, e3, b1), null,
			1, 2, 1, B("GTAC"), new byte[] {1,2,3,4});
		AssemblyAttributes e = new AssemblyAttributes(r);
		assertTrue(e.getOriginatingFragmentID(Range.closed(0,0), ImmutableSet.of(0), null).contains("r1"));
		assertTrue(e.getOriginatingFragmentID(Range.closed(0,0), ImmutableSet.of(0), null).contains("r2"));
		assertTrue(e.getOriginatingFragmentID(Range.closed(0,0), ImmutableSet.of(0), null).contains("r3"));
		assertTrue(e.getOriginatingFragmentID(Range.closed(0,0), ImmutableSet.of(1), null).contains("r4"));

		assertFalse(e.getOriginatingFragmentID(Range.closed(0,0), ImmutableSet.of(1), null).contains("r1"));
		assertFalse(e.getOriginatingFragmentID(Range.closed(0,0), ImmutableSet.of(1), null).contains("r2"));
		assertFalse(e.getOriginatingFragmentID(Range.closed(0,0), ImmutableSet.of(1), null).contains("r3"));
		assertFalse(e.getOriginatingFragmentID(Range.closed(0,0), ImmutableSet.of(0), null).contains("r4"));
	}
	@Test
	public void getEvidenceIDs_should_return_underlying_evidence() {
		DirectedEvidence e1 = SCE(BWD, Read(0, 1, "5S5M"));
		DirectedEvidence e2 = SCE(BWD, Read(0, 1, "6S5M"));
		DirectedEvidence e3 = NRRP(OEA(0, 1, "1M", false));
		SAMRecord e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), FWD, Lists.newArrayList(e1, e2, e3), null,
			1, 2, 1, B("GTAC"), new byte[] {1,2,3,4});
		//Collection<String> ids = e.getEvidenceIDs();
		Collection<String> eid = new AssemblyAttributes(e).getEvidenceIDs(Range.closed(0,0), null, null);
		assertEquals(3, eid.size());
		assertTrue(eid.contains(e1.getEvidenceID()));
		assertTrue(eid.contains(e2.getEvidenceID()));
		assertTrue(eid.contains(e3.getEvidenceID()));
	}
	@Test
	public void getEvidenceIDs_should_return_empty_collection_for_no_evidence() {
		SAMRecord e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), FWD, Lists.newArrayList(),null,
			1, 2, 1, B("GTAC"), new byte[] {1,2,3,4});
		assertEquals(0, new AssemblyAttributes(e).getEvidenceIDs(Range.closed(0, 4), null, null).size());
	}
	@Rule
    public TemporaryFolder testFolder = new TemporaryFolder();
	@Test
	public void realign_should_fix_778_chr1_170849702_misalignment() throws ConfigurationException {
		String assembly = "ATCCATCCCTATGACCCAAACATCTCCCACCAGGCCTCATGTTCAATATTAAAGATCACATTTCAACTTGAGATTTGGAGGGGACAAACATACAAATCATATCATTATCTCTCTCCCCACTTCTCTCTTTATCAATCCCTCCCTCTTTGTCAATCTTAGCCTTGGCCTTCAGATTTTACCACTTGATTTTTCACATTTTCTGTATTCTTAAT"
				+ "GATTATTATATTTTCATGTTCTTGCTAATCTATATCATGGTTAGAAATCAAAGCATGCCGAAATTTCTCTCTTACTTTTTTTGCTGTT";
		File ref = new File("src/test/resources/chr1_170849600_170849850.fa");
		ProcessingContext context = new ProcessingContext(getFSContext(), ref, null, new ArrayList<Header>(),
				new GridssConfiguration());
		context.registerCategory("input.bam");
		SAMRecord e = AssemblyFactory.createAnchoredBreakend(context, AES(context), new SequentialIdGenerator("asm"), BWD, null, null,
				0, 170849702-170849600+1, 97, B(assembly), B(40, assembly.length()));
		e = SAMRecordUtil.realign(context.getReference(), e, 50, true);
		// anchor location is 11bp off
		assertEquals("212S88M", e.getCigarString());
		assertEquals(170849713-170849600+1, asAssemblyEvidence(e).getBreakendSummary().start);
	}
	@Test
	public void small_indel_should_be_called_if_realignment_spans_event() {
		String assembly = "AAAAAAAAAATTAAAAAAAAAA";
		SAMRecord e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), FWD, null,null,
				0, 10, 10, B(assembly), B(40, assembly.length()));
		e = SAMRecordUtil.realign(getContext().getReference(), e, 50, true);
		assertEquals("10M2I10M", e.getCigarString());
	}
	/**
	 * Needed for realignment that turn out to match reference
	 */
	@Test
	public void should_allow_reference_allele_assemblies() {
		String assembly = "AAAAAAAAAAA";
		SAMRecord e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), FWD, null,null,
				0, 1, assembly.length(), B(assembly), B(40, assembly.length()));
		
		assertEquals(0, SingleReadEvidence.createEvidence(SES(), 0, e).size());
	}
	@Test
	public void should_allow_realignment_to_reference_allele() {
		String assembly = "AAAAAAAAAAA";
		SAMRecord e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), FWD, null,null,
				0, 1, 1, B(assembly), B(40, assembly.length()));
		e = SAMRecordUtil.realign(getContext().getReference(), e, 50, true);
		assertEquals(0, SingleReadEvidence.createEvidence(SES(), 0, e).size());
	}
	@Test
	public void getAnchor_should_not_include_inexact_breakend_bases() {
		//List<DirectedEvidence> support = Lists.<DirectedEvidence>newArrayList(NRRP(OEA(0, 1, "1M", true)));
		assertEquals(0, asAssemblyEvidence(AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), new BreakendSummary(0, FWD, 1), null, null, B("GTAC"), new byte[] {1,2,3,4})).getAnchorSequence().length);
	}
	@Test
	public void breakpoint_should_use_indel_cigar() {
		//List<DirectedEvidence> support = Lists.<DirectedEvidence>newArrayList();
		// 1234567890   1234567890
		//         MMIIIDDDDMMMM
		//         NNAAA    TTTT
		SAMRecord r = AssemblyFactory.createAnchoredBreakpoint(getContext(), AES(), new SequentialIdGenerator("asm"), null,null,
				0, 10, 2,
				0, 15, 4,
				B("NNAAATTTT"),
				B("ABCDEFGHI"));
		assertEquals(0, (int)r.getReferenceIndex());
		assertEquals(9, r.getAlignmentStart());
		assertEquals("2M3I4D4M", r.getCigarString());
		assertEquals("ABCDEFGHI", S(r.getBaseQualities()));
		assertEquals("NNAAATTTT", S(r.getReadBases()));
	}
	@Test
	public void realign_should_shift_breakend_to_match_reference() {
		SAMRecord e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), FWD, null,null,
				0, 5, 5, B("AAAATTTT"), new byte[] {1,2,3,4,1,2,3,4});
		e = SAMRecordUtil.realign(getContext().getReference(), e, 50, true);
		assertEquals("TTTT", S(asAssemblyEvidence(e).getBreakendSequence()));
		assertEquals(new BreakendSummary(0, FWD, 4), asAssemblyEvidence(e).getBreakendSummary());
		assertEquals("4M4S", e.getCigarString());
	}
	@Test
	public void realign_should_align_to_reference_with_50bp_margin_around_expected_anchor_interval() {
		int margin = 50;
		for (int startpos = 300 - margin; startpos <= 300 + margin; startpos++) {
			String seq = S("N", 50) + S(Arrays.copyOfRange(RANDOM, 299, 399)); // genomic positions 300-400
			SAMRecord e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), BWD, null,null,
					2, startpos, 100, B(seq), B(40, seq.length()));
			e = SAMRecordUtil.realign(getContext().getReference(), e, 50, true);
			assertEquals(300, asAssemblyEvidence(e).getBreakendSummary().start);
			assertEquals(50, asAssemblyEvidence(e).getBreakendSequence().length);
		}
		// FWD breakend
		for (int startpos = 300 - margin; startpos <= 300 + margin; startpos++) {
			String seq = S(Arrays.copyOfRange(RANDOM, 299, 399)) + S("N", 50);
			SAMRecord e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), FWD, null,null,
					2, startpos + 100 - 1, 100, B(seq), B(40, seq.length()));
			e = SAMRecordUtil.realign(getContext().getReference(), e, 50, true);
			assertEquals(399, asAssemblyEvidence(e).getBreakendSummary().start);
			assertEquals(50, asAssemblyEvidence(e).getBreakendSequence().length);
		}
	}
	@Test
	public void realign_should_expand_window_by_breakend_length_to_allow_for_mapping_over_small_indels() {
		int indelSize = 20;
		String seq = "N" + S(Arrays.copyOfRange(RANDOM, 299-indelSize-100, 299-indelSize)) + S(Arrays.copyOfRange(RANDOM, 299, 399)); // genomic positions 300-400
		SAMRecord e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), BWD, null,null,
				2, 300, 100, B(seq), B(40, seq.length()));
		e = SAMRecordUtil.realign(getContext().getReference(), e, 50, true);
		assertEquals("1S100M20D100M", e.getCigarString());
		
		seq = S(Arrays.copyOfRange(RANDOM, 299, 399)) + S(Arrays.copyOfRange(RANDOM, 399+indelSize, 399+indelSize+100)) + "N";
		e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), FWD, null,null,
				2, 399, 100, B(seq), B(40, seq.length()));
		e = SAMRecordUtil.realign(getContext().getReference(), e, 50, true);
		assertEquals("100M20D100M1S", e.getCigarString());
	}
	@Test
	public void realign_should_allow_small_anchor_deletion() {
		String seq = S(B('N', 100)) + S(Arrays.copyOfRange(RANDOM, 0, 100)) + S(Arrays.copyOfRange(RANDOM, 110, 210));
		SAMRecord e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), BWD, null,null,
				2, 1, 210, B(seq), B(40, seq.length()));
		e = SAMRecordUtil.realign(getContext().getReference(), e, 50, true);
		assertEquals("100S100M10D100M", e.getCigarString());
	}
	@Test
	public void realign_should_allow_small_anchor_insertion() {
		String seq = S(B('N', 100)) + S(Arrays.copyOfRange(RANDOM, 0, 100)) + "NNNNNNNNNN" + S(Arrays.copyOfRange(RANDOM, 100, 200));
		SAMRecord e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), BWD, null,null,
				2, 1, 200, B(seq), B(40, seq.length()));
		e = SAMRecordUtil.realign(getContext().getReference(), e, 50, true);
		assertEquals("100S100M10I100M", e.getCigarString());
	}
	@Test
	@Ignore("Not part of SAMRecordUtil realign functionality")
	public void realign_should_abort_if_anchor_turns_into_soft_clip() {
		String seq = S(Arrays.copyOfRange(RANDOM, 0, 10)) + S(Arrays.copyOfRange(RANDOM, 30, 70));
		SAMRecord e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), FWD, null,null,
				2, 1, 10, B(seq), B(40, seq.length()));
		assertEquals("10M40S", e.getCigarString());
		assertEquals("10M40S", SAMRecordUtil.realign(getContext().getReference(), e, 50, true).getCigarString());
	}
	@Test
	@Ignore("Not part of SAMRecordUtil realign functionality")
	public void realign_should_not_flip_tandem_duplication() {
		ProcessingContext pc = getContext();
		String seq = S(Arrays.copyOfRange(RANDOM, 150, 200)) + S(Arrays.copyOfRange(RANDOM, 100, 200)) + S(Arrays.copyOfRange(RANDOM, 100, 150));
		SAMRecord e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), BWD, null,null,
				2, 100, 50, B(seq), B(40, seq.length()));
		
		// breakend alignment is actually better
		SAMRecord e1 = SAMRecordUtil.realign(pc.getReference(), e, 200, true);
		assertEquals("48S102M50S", e1.getCigarString()); // 50S100M50S + homology
		
		// but if we constrain to require anchor bases, we shouldn't realign
		//SAMRecord e2 = SAMRecordUtil.realign(getContext().getReference(), e, 200, true, 0.5f); // 0.5f;
		// assertEquals("150S50M", e2.getCigarString());
	}
	@Test
	public void realign_should_turn_reference_bubble_into_reference_assembly() {
		SAMRecord ass = AssemblyFactory.createAnchoredBreakpoint(getContext(), AES(), new SequentialIdGenerator("asm"), null,null,
				0, 10, 1,
				0, 17, 1,
				B("AAAAAAAA"),
				B("AAAAAAAA"));
		ass = SAMRecordUtil.realign(getContext().getReference(), ass, 50, true);
		assertEquals(0, SingleReadEvidence.createEvidence(SES(), 0, ass).size());
	}
}

