package au.edu.wehi.idsv.debruijn.positional;

import au.edu.wehi.idsv.*;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.sam.SamTags;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Lists;
import com.google.common.collect.Ordering;
import htsjdk.samtools.SAMRecord;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

import static org.junit.Assert.*;


public class NonReferenceContigAssemblerTest extends TestHelper {
	private NonReferenceContigAssembler caller;
	private EvidenceTracker tracker;
	public List<SAMRecord> go(ProcessingContext pc, boolean collapse, DirectedEvidence... input) {
		caller = create(pc, collapse, input);
		List<SAMRecord> assemblies = Lists.newArrayList(caller);
		return assemblies;
	}
	public List<SAMRecord> go(ProcessingContext pc, boolean collapse, BreakendDirection direction, DirectedEvidence... input) {
		caller = create(pc, collapse, direction, input);
		List<SAMRecord> assemblies = Lists.newArrayList(caller);
		return assemblies;
	}
	private NonReferenceContigAssembler create(ProcessingContext pc, boolean collapse, DirectedEvidence... input) {
		return create(pc, collapse, BreakendDirection.Forward, input);
	}
	private NonReferenceContigAssembler create(ProcessingContext pc, boolean collapse, BreakendDirection direction, DirectedEvidence... input) {
		Arrays.sort(input, DirectedEvidence.ByStartEnd);
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, Arrays.stream(input).map(e -> (SAMEvidenceSource)e.getEvidenceSource()).collect(Collectors.toList()), new File("NonReferenceContigAssemblerTest.bam"));
		int maxEvidenceWidth = aes.getMaxConcordantFragmentSize() - aes.getMinConcordantFragmentSize() + 1;
		int maxReadLength = maxReadLength(input);
		int k = pc.getAssemblyParameters().k;
		int maxPathLength = pc.getAssemblyParameters().positional.maxPathLengthInBases(maxReadLength);
		int maxPathCollapseLength = pc.getAssemblyParameters().errorCorrection.maxPathCollapseLengthInBases(maxReadLength);
		tracker = new EvidenceTracker();
		SupportNodeIterator supportIt = new SupportNodeIterator(k, Arrays.stream(input).iterator(), aes.getMaxConcordantFragmentSize(), tracker, pc.getAssemblyParameters().includePairAnchors, pc.getAssemblyParameters().pairAnchorMismatchIgnoreEndBases);
		AggregateNodeIterator agIt = new AggregateNodeIterator(supportIt);
		Iterator<KmerPathNode> pnIt = new PathNodeIterator(agIt, maxPathLength, k);
		if (collapse) {
			pnIt = new PathCollapseIterator(pnIt, k, maxPathCollapseLength, pc.getAssemblyParameters().errorCorrection.maxBaseMismatchForCollapse, pc.getAssemblyParameters().errorCorrection.collapseBubblesOnly, 0);
			pnIt = new PathSimplificationIterator(pnIt, maxPathLength, maxEvidenceWidth);
		}
		caller = new NonReferenceContigAssembler(pnIt, 0, maxEvidenceWidth + maxReadLength + 2, maxReadLength, k, aes, new SequentialIdGenerator("asm"), tracker, "test", direction, null, null);
		return caller;
	}
	@Test
	public void should_call_simple_fwd_SC() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		SoftClipEvidence sce = SCE(FWD, withSequence("ACGTGGTCGACC", Read(0, 5, "6M6S")));
		List<SAMRecord> output = go(pc, true, sce);
		assertEquals(1, output.size());
		assertEquals("ACGTGGTCGACC", S(output.get(0).getReadBases()));
		assertEquals("6M6S", output.get(0).getCigarString());
		assertEquals(0, (int)output.get(0).getReferenceIndex());
		assertEquals(5, output.get(0).getAlignmentStart());
	}
	@Test
	public void should_set_assembly_attributes() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		SoftClipEvidence sce = SCE(FWD, withSequence("ACGTGGTCGACC", Read(0, 5, "6M6S")));
		List<SAMRecord> output = go(pc, true, sce);
		AssemblyAttributes attr = new AssemblyAttributes(output.get(0));
		assertNotNull(attr);
		assertEquals(1, attr.getSupportingReadCount(5, null, null, null));
	}
	@Test
	public void should_call_simple_bwd_SC() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		SoftClipEvidence sce = SCE(BWD, withSequence("ACGTGGTCGACC", Read(0, 10, "6S6M")));
		List<SAMRecord> output = go(pc, true, sce);
		
		assertEquals(1, output.size());
		assertEquals("ACGTGGTCGACC", S(output.get(0).getReadBases()));
		assertEquals("6S6M", output.get(0).getCigarString());
		assertEquals(0, (int)output.get(0).getReferenceIndex());
		assertEquals(10, output.get(0).getAlignmentStart());
	}
	@Test
	public void should_call_overlapping_SC() {
		// 12345678901234567890
		// ACGTGGCTCGACC
		// ------
		//        ------
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		SoftClipEvidence fwd = SCE(BWD, withSequence("ACGTGGCTCGACC", Read(0, 8, "7S6M")));
		SoftClipEvidence bwd = SCE(FWD, withSequence("ACGTGGCTCGACC", Read(0, 1, "6M7S")));
		List<SAMRecord> output = go(pc, true, fwd, bwd);
		assertEquals(1, output.size());
		assertEquals("6M1I1D6M", output.get(0).getCigarString());
	}
	@Test
	public void should_call_simple_fwd_RP() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		SAMEvidenceSource ses = SES(100, 200);
		DiscordantReadPair e = (DiscordantReadPair)NRRP(ses, withSequence("ACGTGGTCGACC", DP(0, 50, "12M", true, 1, 1, "12M", false)));
		List<SAMRecord> output = go(pc, true, e);
		assertEquals(1, output.size());
		SoftClipEvidence ass = SoftClipEvidence.create(AES(ses), FWD, output.get(0));
		assertEquals(e.getBreakendSummary().localBreakend(), ass.getBreakendSummary());
		assertEquals("ACGTGGTCGACC", S(ass.getBreakendSequence()));
		assertEquals("", S(ass.getAnchorSequence()));
	}
	@Test
	public void should_call_simple_bwd_RP() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		MockSAMEvidenceSource ses = SES(100, 200);
		DiscordantReadPair e = (DiscordantReadPair)NRRP(ses, withSequence("ACGTGGTCGACC", DP(0, 450, "12M", false, 1, 1, "12M", true)));
		List<SAMRecord> output = go(pc, true, e);
		assertEquals(1, output.size());
		SoftClipEvidence ass = SoftClipEvidence.create(AES(ses), BWD, output.get(0));
		assertEquals(e.getBreakendSummary().localBreakend(), ass.getBreakendSummary());
		assertEquals("ACGTGGTCGACC", S(ass.getBreakendSequence()));
		assertEquals("", S(ass.getAnchorSequence()));
	}
	@Test
	public void should_handle_sc_repeated_kmers() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		SoftClipEvidence sce = SCE(FWD, withSequence("TTTTTTTTTTTT", Read(0, 5, "6M6S")));
		List<SAMRecord> output = go(pc, true, sce);
		assertEquals(1, output.size());
		assertEquals("TTTTTTTTTTTT", S(output.get(0).getReadBases()));
	}
	@Test
	public void should_call_multiple_nonoverlapping_contigs() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		SoftClipEvidence sce = SCE(FWD, withSequence("ACGTGGTCGACC", Read(0, 5, "6M6S")));
		DiscordantReadPair e = (DiscordantReadPair)NRRP(SES(100, 200), withSequence("GACCTCTACT", DP(0, 25, "10M", true, 1, 1, "10M", false)));
		List<SAMRecord> output = go(pc, true, sce, e);
		assertEquals(2, output.size());
		assertEquals("ACGTGGTCGACC", S(output.get(0).getReadBases()));
		assertEquals("NN" + "GACCTCTACT", S(output.get(1).getReadBases()));
	}
	@Test
	public void should_call_overlapping_sc_rp() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		SoftClipEvidence sce = SCE(FWD, withSequence("ACGTGGTCGACC", Read(0, 50, "6M6S")));
		DiscordantReadPair e = (DiscordantReadPair)NRRP(SES(10, 200), withSequence("GACCTCTACT", DP(0, 25, "10M", true, 1, 1, "10M", false)));
		List<SAMRecord> output = go(pc, true, sce, e);
		assertEquals(1, output.size());
		assertEquals("ACGTGGTCGACCTCTACT", S(output.get(0).getReadBases()));
	}
	@Test
	public void should_remove_fully_reference_evidence() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		SoftClipEvidence sce =     SCE(FWD, withSequence("ACGTGGTCGACC", Read(0, 50, "6M6S")));
		SoftClipEvidence fullref = SCE(FWD, withSequence("ACGTGG", Read(0, 50, "5M1S")));
		List<SAMRecord> output = go(pc, true, sce, fullref);
		assertEquals(1, output.size());
		assertEquals("ACGTGGTCGACC", S(output.get(0).getReadBases()));
		assertEquals(2, tracker.tracking_evidenceTotal());
		assertEquals(0, tracker.tracking_evidenceActive());
		assertEquals(0, caller.tracking_activeNodes());
	}
	@Test
	public void should_remove_fully_reference_evidence_before_end() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		int n = 8;
		DirectedEvidence[] input = new DirectedEvidence[2 * n];
		for (int i = 0; i < n; i++) {
			input[2*i] =      SCE(FWD, withSequence( "ACGTGGTCGACC", Read(0, 1000*i, "6M6S")));
			input[2*i + 1 ] = SCE(FWD, withSequence("AACGTGG", Read(0, 1000*i-1, "6M1S")));
		}
		NonReferenceContigAssembler caller = create(pc, true, input);
		for (int i = 0; i < n / 2; i++) {
			caller.next();
		}
		assertTrue(caller.tracking_activeNodes() <= 3);
	}
	@Test
	public void should_truncate_anchored_polyA_expansion_fwd() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		SoftClipEvidence sce = SCE(FWD, withSequence("GTACAAAA", Read(0, 10, "4M4S")));
		DiscordantReadPair dp = (DiscordantReadPair)NRRP(SES(11, 100), withSequence("CAAAAT", DP(0, 1, "6M", true, 1, 1, "6M", false)));
		List<SAMRecord> output = go(pc, true, sce, dp);
		assertEquals(1, output.size());
		assertEquals("GTACAAAA", S(output.get(0).getReadBases()));
	}
	@Test
	public void should_truncate_anchored_polyA_expansion_bwd() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		SoftClipEvidence sce = SCE(BWD, withSequence("AAAAGTAC", Read(0, 100, "4S4M")));
		DiscordantReadPair dp = (DiscordantReadPair)NRRP(SES(11, 100), withSequence("CAAAAG", DP(0, 110, "6M", false, 1, 1, "6M", true)));
		List<SAMRecord> output = go(pc, true, sce, dp);
		assertEquals(1, output.size());
		assertEquals("AAAAGTAC", S(output.get(0).getReadBases()));
	}
	@Test
	public void should_take_higher_weight_unanchored_assembly() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		pc.getAssemblyParameters().includePairAnchors = false;
		List<SAMRecord> output = go(pc, true,
				NRRP(SES(1, 100), withSequence("CAAAAAAAAAAAGT", DP(1, 100, "14M", false, 1, 1, "14M", true))),
				NRRP(SES(1, 100), withSequence("CAAAAAAAAAAGTC", DP(1, 101, "14M", false, 1, 1, "14M", true))));
		assertEquals(1, output.size());
		assertEquals("AAAAAAAAAAAGTC" + "NN", S(output.get(0).getReadBases()));
	}
	@Test
	public void should_handle_ambigious_bases() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		pc.getAssemblyParameters().includePairAnchors = false;
		List<SAMRecord> output = go(pc, true,
				NRRP(SES(1, 100), withSequence("AATAAGAGT", DP(1, 100, "9M", false, 1, 1, "9M", true))),
				NRRP(SES(1, 100), withSequence("AATGNGAGTC", DP(1, 101, "10M", false, 1, 1, "10M", true))));
		assertEquals(1, output.size());
		assertEquals("AATAAGAGTC" + "NN", S(output.get(0).getReadBases()));
	}
	@Test
	public void should_handle_multiple_candidate_offsets_for_a_single_kmer_position() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		pc.getAssemblyParameters().includePairAnchors = false;
		List<SAMRecord> output = go(pc, true,
				NRRP(SES(1, 100), withSequence("TAAAAA", DP(1, 100, "6M", false, 0, 1, "6M", true))));
		assertEquals(1, output.size());
		assertEquals("TAAAAA" + "NN", S(output.get(0).getReadBases()));
	}
	@Test
	public void fwd_anchor_length_calculation_should_not_include_any_nonreference_bases() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 2;
		List<SAMRecord> output = go(pc, true,
				SCE(FWD, withSequence("AAAAAAAAAT", Read(0, 1, "9M1S"))),
				SCE(FWD, withSequence(        "ATTTTTTTTT", Read(0, 9, "1M9S"))));
		assertEquals(1, output.size());
		assertEquals("AAAAAAAAATTTTTTTTT", S(output.get(0).getReadBases()));
	}
	@Test
	public void bwd_anchor_length_calculation_should_not_include_any_nonreference_bases() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 2;
		List<SAMRecord> output = go(pc, true,
				SCE(BWD, withSequence("AAAAAAAAAT", Read(0, 10, "9S1M"))),
				SCE(BWD, withSequence(        "ATTTTTTTTT", Read(0, 10, "1S9M"))));
		assertEquals(1, output.size());
		assertEquals("AAAAAAAAATTTTTTTTT", S(output.get(0).getReadBases()));
	}
	@Test
	public void anchor_should_be_trimmed_at_max_length() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 3;
		List<SAMRecord> output = go(pc, true,
				SCE(FWD, withSequence("AAAAAAAAAT", Read(0, 1, "9M1S"))),
				SCE(FWD, withSequence("AAAAAAAAAT", Read(0, 2, "9M1S"))),
				SCE(FWD, withSequence("AAAAAAAAAT", Read(0, 3, "9M1S"))),
				SCE(FWD, withSequence("AAAAAAAAAT", Read(0, 4, "9M1S"))),
				SCE(FWD, withSequence("AAAAAAAAAT", Read(0, 5, "9M1S"))),
				SCE(FWD, withSequence("AAAAAAAAAT", Read(0, 6, "9M1S"))),
				SCE(FWD, withSequence("AAAAAAAAAT", Read(0, 7, "9M1S"))),
				SCE(FWD, withSequence("AAAAAAAAAT", Read(0, 8, "9M1S"))),
				SCE(FWD, withSequence("AAAAAAAAAT", Read(0, 9, "9M1S"))),
				SCE(FWD, withSequence("AAAAAAAAAT", Read(0, 10, "9M1S"))),
				SCE(FWD, withSequence("AAAAAAAAAT", Read(0, 11, "9M1S"))),
				SCE(FWD, withSequence("AAAAAAAAAT", Read(0, 12, "9M1S"))),
				SCE(FWD, withSequence("AAAAAAAAAT", Read(0, 13, "9M1S"))),
				SCE(FWD, withSequence("AAAAAAAAAT", Read(0, 14, "9M1S"))),
				SCE(FWD, withSequence(       "AATTTTTTTTT", Read(0, 18, "2M9S"))));
		assertEquals("AAAAAAAAAAATTTTTTTTT", S(output.get(0).getReadBases()));
	}
	@Test
	public void should_assembly_UM_into_anchor() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		pc.getAssemblyParameters().includePairAnchors = true;
		pc.getAssemblyParameters().pairAnchorMismatchIgnoreEndBases = 0;
		SAMRecord[] rp = OEA(0, 100, "10M", true);
		rp[0].setReadBases(B(S(RANDOM).substring(0, 10)));
		rp[1].setReadBases(B(S(RANDOM).substring(2, 12)));
		rp[1].setReadNegativeStrandFlag(true);
		NonReferenceReadPair nrrp = NRRP(SES(1, 100), rp);
		List<SAMRecord> output = go(pc, false, nrrp);
		assertEquals(S(RANDOM).substring(0, 12), S(output.get(0).getReadBases()));
		assertEquals(0, (int)output.get(0).getReferenceIndex());
		assertEquals(100, output.get(0).getAlignmentStart());
		assertEquals("10M2S", output.get(0).getCigarString());
	}
	@Test
	public void should_removeMisassembledPartialContigsDuringAssembly() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().removeMisassembledPartialContigsDuringAssembly = false;
		MockSAMEvidenceSource ses = new MockSAMEvidenceSource(getContext(), 50, 100);
		pc.getAssemblyParameters().k = 4;
		pc.getAssemblyParameters().maxExpectedBreakendLengthMultiple = 1.5f;
		List<DirectedEvidence> e = new ArrayList<>();
		for (int i = 1; i < 1000; i++) {
			e.add(NRRP(ses, OEA(0, i, "10M", true)));
		}
		List<SAMRecord> output = go(pc, false, e.toArray(new DirectedEvidence[0]));
		assertEquals(1, output.size());
		SAMRecord r = output.get(0);
		int sclenFull = SAMRecordUtil.getStartClipLength(r) + SAMRecordUtil.getEndClipLength(r);
		
		pc.getAssemblyParameters().removeMisassembledPartialContigsDuringAssembly = true;
		
		output = go(pc, false, e.toArray(new DirectedEvidence[0]));
		assertEquals(1, output.size());
		r = output.get(0);
		int sclenTruncated = SAMRecordUtil.getStartClipLength(r) + SAMRecordUtil.getEndClipLength(r);
		Assert.assertNotEquals(sclenFull, sclenTruncated);
	}
	@Test
	public void should_call_suboptimal_contigs_if_graph_size_limit_exceeded() {
		ProcessingContext pc = getContext();
		MockSAMEvidenceSource ses = SES(10, 10);
		pc.getAssemblyParameters().k = 4;
		pc.getAssemblyParameters().maxExpectedBreakendLengthMultiple = 1;
		pc.getAssemblyParameters().positional.retainWidthMultiple = 1;
		pc.getAssemblyParameters().positional.flushWidthMultiple = 2;
		List<DirectedEvidence> e = new ArrayList<>();
		for (int i = 0; i < 100; i++) {
			for (int j = 0; j < i + 1; j++) {
				e.add(SCE(FWD, ses, withReadName(String.format("%d-%d", i, j), withSequence("AACCGGTT", Read(0, i, "4M4S")))[0]));
			}
		}
		List<SAMRecord> output = go(pc, false, e.toArray(new DirectedEvidence[0]));
		assertEquals(100, output.size());
		// unbounded execution would call high->low thus call the last contig first
		// bounded execution will force the early calls to be made first
		List<Integer> startPos = output.stream().map(x -> x.getAlignmentStart()).collect(Collectors.toList());
		assertFalse(Ordering.natural().reverse().isOrdered(startPos));
	}
	@Test
	public void should_not_assembly_when_graph_density_exceeds_maximum() {
		ProcessingContext pc = getContext();
		MockSAMEvidenceSource ses = SES(10, 10);
		pc.getAssemblyParameters().k = 4;
		pc.getAssemblyParameters().maxExpectedBreakendLengthMultiple = 1;
		pc.getAssemblyParameters().positional.maximumNodeDensity = 0.1f;
		List<DirectedEvidence> e = new ArrayList<>();
		for (int i = 1; i < 101; i++) {
			e.add(SCE(FWD, ses, withReadName(String.format("%d-%d", i, 0), withSequence("AAAATTGG", Read(0, i, "4M4S")))[0]));
			e.add(SCE(FWD, ses, withReadName(String.format("%d-%d", i, 1), withSequence("AAAACCGG", Read(0, i, "4M4S")))[0]));
		}
		List<SAMRecord> output = go(pc, false, e.toArray(new DirectedEvidence[0]));
		// the final advance doesn't get flushed since the underlying stream
		// is considered to advance to Integer.MAX_VALUE at end of input
		assertTrue(output.size() <= 20);
		
		// shouldn't get throttled
		pc.getAssemblyParameters().positional.maximumNodeDensity = 10;
		output = go(pc, false, e.toArray(new DirectedEvidence[0]));
		assertEquals(2 * 100, output.size());
	}
	@Test
	public void should_remove_misassembled_partial_paths() {
		ProcessingContext pc = getContext();
		MockSAMEvidenceSource ses = SES(10, 10);
		pc.getAssemblyParameters().k = 4;
		pc.getAssemblyParameters().maxExpectedBreakendLengthMultiple = 1;
		List<DirectedEvidence> e = new ArrayList<>();
		for (int i = 1; i < 100; i++) {
			e.add(SCE(FWD, ses, withReadName(String.format("%d-%d", 0, i), withSequence("TAAAAAAAAAAAAAAAAAAAAAAAAAAAA", Read(0, i*20, "4M25S")))[0]));
		}
		List<SAMRecord> output = go(pc, false, e.toArray(new DirectedEvidence[0]));
		// the rest should have been removed
		assertEquals(1, output.size());
	}
	@Test
	public void should_not_remove_misassembled_partial_paths_when_could_be_collapsed_with_adjacent_read() {
		ProcessingContext pc = getContext();
		MockSAMEvidenceSource ses = SES(10, 10);
		pc.getAssemblyParameters().k = 2;
		pc.getAssemblyParameters().maxExpectedBreakendLengthMultiple = 5;
		pc.getAssemblyParameters().removeMisassembledPartialContigsDuringAssembly = true;
		pc.getAssemblyParameters().positional.flushWidthMultiple = 100000;
		pc.getAssemblyParameters().positional.maxPathLengthMultiple = 4;
		List<DirectedEvidence> e = new ArrayList<>();
		for (int i = 1; i < 1000; i += 1) {
			SAMRecord[] rp = OEA(0, i, "2M", true);
			rp[0].setReadBases(B("GG"));
			rp[1].setReadBases(B("AA")); // expect reverse comp of read sequence
			e.add(NRRP(ses, rp));
		}
		go(pc, true, e.toArray(new DirectedEvidence[0]));
	}
	@Test
	public void should_generate_support_attributes() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		MockSAMEvidenceSource ses0 = SES(false);
		MockSAMEvidenceSource ses1 = SES(true);
		SoftClipEvidence sce1 = SCE(FWD, ses0, withSequence("ACGTGGTCGACC", Read(0, 5, "6M6S")));
		SoftClipEvidence sce2 = SCE(FWD, ses0, withSequence("ACGTGGTCGAC", Read(0, 5, "6M5S")));
		SoftClipEvidence sce3 = SCE(FWD, ses1, withSequence("ACGTGGTCGT", Read(0, 5, "6M4S")));
		List<SAMRecord> output = go(pc, true, sce1, sce2, sce3);
		assertEquals(1, output.size());
		assertEquals(3, output.get(0).getSignedByteArrayAttribute(SamTags.ASSEMBLY_EVIDENCE_TYPE).length);
		assertEquals(3, output.get(0).getSignedIntArrayAttribute(SamTags.ASSEMBLY_EVIDENCE_CATEGORY).length);
		assertEquals(3, output.get(0).getSignedIntArrayAttribute(SamTags.ASSEMBLY_EVIDENCE_OFFSET_START).length);
		assertEquals(3, output.get(0).getSignedIntArrayAttribute(SamTags.ASSEMBLY_EVIDENCE_OFFSET_END).length);
		assertEquals(3, output.get(0).getFloatArrayAttribute(SamTags.ASSEMBLY_EVIDENCE_QUAL).length);
		assertNotNull(output.get(0).getAttribute(SamTags.ASSEMBLY_EVIDENCE_EVIDENCEID));
		assertNotNull(output.get(0).getAttribute(SamTags.ASSEMBLY_EVIDENCE_FRAGMENTID));

	}
	@Test
	public void should_generate_rp_contig_support_over_rp_interval() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 5;
		pc.getAssemblyParameters().pairAnchorMismatchIgnoreEndBases = 0;
		MockSAMEvidenceSource ses = new MockSAMEvidenceSource(pc, 0, 30);
		
		// 123456789012345678901234567
		// GGTTGCATAGACGTGGTCGACCTAGTA
		// ==========       ==========
		//      MMMMMSSSSSSSSSSSS
		SAMRecord[] rp = DP(0, 1, "10M", true, 0, 1000, "10M", false);
		rp[0].setReadBases(B("GGTTGCATAG"));
		rp[1].setReadBases(B("CGACCTAGTA"));
		SoftClipEvidence sce = SCE(FWD, ses, withSequence("CATAGACGTGGTCGACC", Read(0, 6, "5M12S")));
		List<SAMRecord> output = go(pc, true, sce, NonReferenceReadPair.create(rp[0], rp[1], ses));
		assertEquals(1, output.size());
		assertEquals(27, output.get(0).getReadLength());
		assertEquals("GGTTGCATAGACGTGGTCGACCTAGTA", S(output.get(0).getReadBases()));
		AssemblyAttributes aa = new AssemblyAttributes(output.get(0));
		// GGTTGCATAGACGTGGTCGACCTAGTA
		// 012345678901234567890123456
		assertEquals(0, aa.getSupportingReadCount(0, null, ImmutableSet.of(AssemblyEvidenceSupport.SupportType.ReadPair), null));
		for (int i = 1; i <= 26; i++) {
			assertEquals(1, aa.getSupportingReadCount(i, null, ImmutableSet.of(AssemblyEvidenceSupport.SupportType.ReadPair), null));
		}
		assertEquals(0, aa.getSupportingReadCount(27, null, ImmutableSet.of(AssemblyEvidenceSupport.SupportType.ReadPair), null));
	}
	@Test
	public void rp_support_without_anchor_should_extend_to_contig_start() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 5;
		pc.getAssemblyParameters().pairAnchorMismatchIgnoreEndBases = 0;
		MockSAMEvidenceSource ses = new MockSAMEvidenceSource(pc, 0, 30);
		SAMRecord[] rp = DP(0, 1, "10M", true, 0, 1000, "10M", false);
		rp[0].setReadBases(B("NNNNNNNNNN"));
		rp[1].setReadBases(B("CGACCTAGTA"));
		SoftClipEvidence sce = SCE(FWD, ses, withSequence("CATAGACGTGGTCGACC", Read(0, 6, "5M12S")));
		List<SAMRecord> output = go(pc, true, sce, NonReferenceReadPair.create(rp[0], rp[1], ses));
		assertEquals(1, output.size());
		assertEquals("CATAGACGTGGTCGACCTAGTA", S(output.get(0).getReadBases()));
		AssemblyAttributes aa = new AssemblyAttributes(output.get(0));
		// CATAGACGTGGTCGACCTAGTA
		// 0123456789012345678901
		for (int i = 0; i <= 21; i++) {
			assertEquals(1, aa.getSupportingReadCount(i, null, ImmutableSet.of(AssemblyEvidenceSupport.SupportType.ReadPair), null));
		}
		assertEquals(0, aa.getSupportingReadCount(22, null, ImmutableSet.of(AssemblyEvidenceSupport.SupportType.ReadPair), null));
	}

	@Test
	public void rp_support_without_anchor_should_extend_to_contig_end() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 5;
		pc.getAssemblyParameters().pairAnchorMismatchIgnoreEndBases = 0;
		MockSAMEvidenceSource ses = new MockSAMEvidenceSource(pc, 0, 100);
		SAMRecord[] rp = DP(0, 120, "10M", false, 0, 1000, "10M", true);
		rp[0].setReadBases(B("NNNNNNNNNN"));
		rp[1].setReadBases(B("GGTTGCATAG"));
		SoftClipEvidence sce = SCE(BWD, ses, withSequence("CATAGACGTGGTCGACC", Read(0, 100, "12S5M")));
		List<SAMRecord> output = go(pc, true, BWD, sce, NonReferenceReadPair.create(rp[0], rp[1], ses));
		assertEquals(1, output.size());
		assertEquals("GGTTGCATAGACGTGGTCGACC", S(output.get(0).getReadBases()));
		assertFalse(output.get(0).hasAttribute(SamTags.UNANCHORED));
		AssemblyAttributes aa = new AssemblyAttributes(output.get(0));
		// GGTTGCATAGACGTGGTCGACC
		// 0123456789012345678901
		assertEquals(0, aa.getSupportingReadCount(0, null, ImmutableSet.of(AssemblyEvidenceSupport.SupportType.ReadPair), null));
		for (int i = 1; i < output.get(0).getReadLength(); i++) {
			assertEquals(1, aa.getSupportingReadCount(i, null, ImmutableSet.of(AssemblyEvidenceSupport.SupportType.ReadPair), null));
		}
		assertEquals(0, aa.getSupportingReadCount(output.get(0).getReadLength(), null, ImmutableSet.of(AssemblyEvidenceSupport.SupportType.ReadPair), null));
	}
	@Test
	public void should_not_write_inconsistent_anchors() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		MockSAMEvidenceSource ses = new MockSAMEvidenceSource(pc, 0, 100);
		SAMRecord[] rp = DP(0, 1, "10M", true, 0, 1000, "10M", false);
		rp[0].setReadBases(B("NNNNNNNNNN"));
		rp[1].setReadBases(B("ACGTTGGCCA"));
		rp[1].setReadUnmappedFlag(true);
		rp[0].setMateUnmappedFlag(true);
		SoftClipEvidence sce1 = SCE(FWD, ses, withSequence("ACGTAAAA", Read(0, 20, "4M4S")));
		SoftClipEvidence sce2 = SCE(FWD, ses, withSequence("ACGTAAAA", Read(0, 21, "4M4S")));
		SoftClipEvidence sce3 = SCE(FWD, ses, withSequence("ACGTAAAA", Read(0, 22, "4M4S")));
		SoftClipEvidence sce4 = SCE(FWD, ses, withSequence("ACGTAAAA", Read(0, 23, "4M4S")));
		SoftClipEvidence sce5 = SCE(FWD, ses, withSequence("GCCAAAAA", Read(0, 40, "4M4S")));
		List<SAMRecord> output = go(pc, true, sce1, sce2, sce3, sce4, sce5, NonReferenceReadPair.create(rp[0], rp[1], ses));
		assertFalse(output.get(0).getCigarString().contains("I"));
		assertFalse(output.get(0).getCigarString().contains("D"));
	}
	@Test
	public void unanchored_rp_should_support_at_contig_bounds_but_not_padding_bases_FWD() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 5;
		pc.getAssemblyParameters().minReads = 1;
		pc.getAssemblyParameters().pairAnchorMismatchIgnoreEndBases = 0;
		MockSAMEvidenceSource ses = new MockSAMEvidenceSource(pc, 0, 50);
		SAMRecord[] rp = DP(0, 1, "10M", true, 0, 1000, "10M", false);
		rp[0].setReadBases(B("NNNNNNNNNN"));
		rp[1].setReadBases(B("CGACCTAGTA"));
		List<SAMRecord> output = go(pc, true, FWD, NonReferenceReadPair.create(rp[0], rp[1], ses));
		assertEquals(1, output.size());
		assertEquals(rp[1].getReadLength() + 2, output.get(0).getReadLength());
		AssemblyAttributes aa = new AssemblyAttributes(output.get(0));
		assertTrue(output.get(0).hasAttribute(SamTags.UNANCHORED));
		assertEquals("1X29N1X10S", output.get(0).getCigarString());
		// 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4
		//  * * M M M M M M M M M M - -
		assertEquals(0, aa.getSupportingReadCount(0, null, null, null));
		assertEquals(0, aa.getSupportingReadCount(1, null, null, null));
		assertEquals(1, aa.getSupportingReadCount(2, null, null, null));
		assertEquals(1, aa.getSupportingReadCount(3, null, null, null));
		assertEquals(1, aa.getSupportingReadCount(4, null, null, null));
		assertEquals(1, aa.getSupportingReadCount(10, null, null, null));
		assertEquals(1, aa.getSupportingReadCount(11, null, null, null));
		assertEquals(0, aa.getSupportingReadCount(12, null, null, null));
		assertEquals(0, aa.getSupportingReadCount(13, null, null, null));
	}
	@Test
	public void unanchored_rp_should_support_at_contig_bounds_but_not_padding_bases_BWD() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 5;
		pc.getAssemblyParameters().minReads = 1;
		pc.getAssemblyParameters().pairAnchorMismatchIgnoreEndBases = 0;
		MockSAMEvidenceSource ses = new MockSAMEvidenceSource(pc, 0, 50);
		SAMRecord[] rp = DP(0, 100, "10M", false, 0, 1000, "10M", false);
		rp[0].setReadBases(B("NNNNNNNNNN"));
		rp[1].setReadBases(B("CGACCTAGTA"));
		List<SAMRecord> output = go(pc, true, BWD, NonReferenceReadPair.create(rp[0], rp[1], ses));
		assertEquals(1, output.size());
		assertEquals(rp[1].getReadLength() + 2, output.get(0).getReadLength());
		AssemblyAttributes aa = new AssemblyAttributes(output.get(0));
		assertTrue(output.get(0).hasAttribute(SamTags.UNANCHORED));
		assertEquals("10S1X29N1X", output.get(0).getCigarString());
		assertEquals(0, aa.getSupportingReadCount(0, null, null, null));
		//  0 1 2 3 4 5 6 7 8 9 0 1 2 3
		// M M M M M M M M M M * * - -
		assertEquals(0, aa.getSupportingReadCount(0, null, null, null));
		assertEquals(1, aa.getSupportingReadCount(1, null, null, null));
		assertEquals(1, aa.getSupportingReadCount(2, null, null, null));
		assertEquals(1, aa.getSupportingReadCount(3, null, null, null));
		assertEquals(1, aa.getSupportingReadCount(8, null, null, null));
		assertEquals(1, aa.getSupportingReadCount(9, null, null, null));
		assertEquals(0, aa.getSupportingReadCount(10, null, null, null));
		assertEquals(0, aa.getSupportingReadCount(11, null, null, null));
		assertEquals(0, aa.getSupportingReadCount(12, null, null, null));
		assertEquals(0, aa.getSupportingReadCount(13, null, null, null));
	}
	@Test
	public void issue287_rp_support_for_anchor_with_no_valid_kmers_treated_as_anchoring_read_not_assembled() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 5;
		pc.getAssemblyParameters().pairAnchorMismatchIgnoreEndBases = 0;
		MockSAMEvidenceSource ses = new MockSAMEvidenceSource(pc, 0, 100);
		SAMRecord[] rp = DP(0, 120, "1M", false, 0, 1000, "10M", true);
		rp[0].setReadBases(B("A"));
		rp[1].setReadBases(B("GGTTGCATAG"));
		SoftClipEvidence sce = SCE(BWD, ses, withSequence("CATAGACGTGGTCGACC", Read(0, 100, "12S5M")));
		List<SAMRecord> output = go(pc, true, BWD, sce, NonReferenceReadPair.create(rp[0], rp[1], ses));
		assertEquals(1, output.size());
		assertEquals("GGTTGCATAGACGTGGTCGACC", S(output.get(0).getReadBases()));
		assertFalse(output.get(0).hasAttribute(SamTags.UNANCHORED));
		AssemblyAttributes aa = new AssemblyAttributes(output.get(0));
		// GGTTGCATAGACGTGGTCGACC
		// 0123456789012345678901
		assertEquals(0, aa.getSupportingReadCount(0, null, ImmutableSet.of(AssemblyEvidenceSupport.SupportType.ReadPair), null));
		for (int i = 1; i < output.get(0).getReadLength(); i++) {
			assertEquals(1, aa.getSupportingReadCount(i, null, ImmutableSet.of(AssemblyEvidenceSupport.SupportType.ReadPair), null));
		}
		assertEquals(0, aa.getSupportingReadCount(output.get(0).getReadLength(), null, ImmutableSet.of(AssemblyEvidenceSupport.SupportType.ReadPair), null));
	}
}
