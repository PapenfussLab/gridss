package au.edu.wehi.idsv.debruijn.positional;

import au.edu.wehi.idsv.*;
import au.edu.wehi.idsv.bed.IntervalBed;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Lists;
import htsjdk.samtools.SAMRecord;
import org.junit.Assert;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.junit.Assert.*;


public class PositionalAssemblerTest extends TestHelper {
	@Test
	public void should_assemble_simple_input() {
		ProcessingContext pc = getContext();
		AssemblyEvidenceSource aes = AES(pc);
		pc.getAssemblyParameters().k = 4;
		List<DirectedEvidence> input = new ArrayList<DirectedEvidence>();
		input.add(SCE(BWD, Read(0, 10, "5S5M")));
		input.add(SCE(FWD, Read(0, 10, "5M5S")));
		input.add(SCE(BWD, Read(0, 100, "5S5M")));
		input.add(SCE(FWD, Read(0, 100, "5M5S")));
		input.sort(DirectedEvidenceOrder.ByStartEndStart2End2);
		ArrayList<SAMRecord> r = Lists.newArrayList(new PositionalAssembler(pc, aes, new SequentialIdGenerator("asm"), input.iterator(), null, null));
		assertEquals(4, r.size());
	}
	@Test
	public void should_assemble_simple_forward_soft_clips() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().anchorLength = 1;
		AssemblyEvidenceSource aes = AES(pc);
		pc.getAssemblyParameters().k = 4;
		List<DirectedEvidence> input = new ArrayList<DirectedEvidence>();
		// 12345678901234567890
		//          ACGTTGGTTA
		//          MMMMMSSSSS
		//       GCAACGTTGGTTAA
		//       MMMMMMMMSSSSSS
		input.add(SCE(FWD, withSequence(        "ACGTTGGTTA", Read(0, 10, "5M5S"))[0]));
		input.add(SCE(FWD, withSequence("TTTTTGCAACGTTGGTTAA", Read(0, 2, "13M6S"))[0]));
		input.sort(DirectedEvidenceOrder.ByStartEndStart2End2);
		List<SingleReadEvidence> r = asAssemblyEvidence(aes, Lists.newArrayList(new PositionalAssembler(pc, aes, new SequentialIdGenerator("asm"), input.iterator(), null, null)));
		assertEquals(1, r.size());
		assertEquals(new BreakendSummary(0, FWD, 14), r.get(0).getBreakendSummary());
		// anchor length to match breakend length
		assertEquals("AACGTT", S(r.get(0).getAnchorSequence()));
		assertEquals("GGTTAA", S(r.get(0).getBreakendSequence()));
		assertEquals("AACGTTGGTTAA", S(r.get(0).getSAMRecord().getReadBases()));
        assertTrue(AssemblyAttributes.isAssembly(r.get(0).getSAMRecord()));
        AssemblyAttributes aa = new AssemblyAttributes(r.get(0).getSAMRecord());
                                  // A A C G T T G G T T A A
        int[] expected = new int[] {1,1,2,2,2,2,2,2,2,2,2,1,0,};
        int[] actual = IntStream.range(0, 12+1).map(i -> aa.getSupportingReadCount(i, null, null, null)).toArray();
		assertArrayEquals(expected, actual);
        assertArrayEquals(expected, IntStream.range(0, 12+1).map(i -> aa.getSupportingReadCount(i, null, ImmutableSet.of(AssemblyEvidenceSupport.SupportType.Read), null)).toArray());
        assertArrayEquals(expected, IntStream.range(0, 12+1).map(i -> aa.getSupportingReadCount(i, ImmutableSet.of(0), ImmutableSet.of(AssemblyEvidenceSupport.SupportType.Read), aes)).toArray());
        assertArrayEquals(new int[13], IntStream.range(0, 12+1).map(i -> aa.getSupportingReadCount(i, ImmutableSet.of(1), null, aes)).toArray());
        assertArrayEquals(new double[13], IntStream.range(0, 12+1).mapToDouble(i -> aa.getSupportingQualScore(i, null, ImmutableSet.of(AssemblyEvidenceSupport.SupportType.ReadPair), null)).toArray(), 0);
	}
	@Test
	public void rp_anchor_should_set_non_reference_bases_as_anchoring() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().anchorLength = 1;
		AssemblyEvidenceSource aes = AES(pc);
		pc.getAssemblyParameters().k = 4;
		List<DirectedEvidence> input = new ArrayList<DirectedEvidence>();
		// 12345678901234567890
		//          ACGTTGGTTA
		//          MMMMMSSSSS
		//       GCAACGTTGGTTAA
		//       MMMMMMMMSSSSSS
		//  TTTTTGCAACGTTGGTTAA
		input.add(SCE(FWD, withSequence("ACGTTGGTTA", Read(0, 10, "5M5S"))[0]));
		input.add(SCE(FWD, withSequence("TTTTTGCAACGTTGGTTAA", Read(0, 2, "13M6S"))[0]));
		input.add(NRRP(withSequence("TTTTTGCAACGTTGGTTAA", DP(0, 2, "13M6S", true, 1, 1, "19M", false))));
		input.sort(DirectedEvidenceOrder.ByStartEndStart2End2);
		List<SingleReadEvidence> r = asAssemblyEvidence(aes, Lists.newArrayList(new PositionalAssembler(pc, aes, new SequentialIdGenerator("asm"), input.iterator(), null, null)));
		assertEquals(2, r.size());
		// race condition w.r.t which assembly returns first
		assertEquals(new BreakendSummary(0, FWD, 14), r.get(0).getBreakendSummary());
		// anchor length to match breakend length
		assertEquals("AACGTT", S(r.get(0).getAnchorSequence()));
		assertEquals("GGTTAA", S(r.get(0).getBreakendSequence()));
		assertEquals("AACGTTGGTTAA", S(r.get(0).getSAMRecord().getReadBases()));
	}
	@Test
	public void anchor_should_not_overrun_contig_start() {
		ProcessingContext pc = getContext();
		AssemblyEvidenceSource aes = AES(pc);
		pc.getAssemblyParameters().k = 4;
		List<DirectedEvidence> input = new ArrayList<DirectedEvidence>();
		input.add(SCE(FWD, withSequence("CGTAACCGGTTC", Read(0, 1, "1M2I4M5S"))));
		ArrayList<SAMRecord> r = Lists.newArrayList(new PositionalAssembler(pc, aes, new SequentialIdGenerator("asm"), input.iterator(), null, null));
		assertEquals(1, r.size());
		assertEquals(1, r.get(0).getAlignmentStart());
	}
	@Test
	public void anchor_should_not_overrun_contig_end() {
		ProcessingContext pc = getContext();
		AssemblyEvidenceSource aes = AES(pc);
		pc.getAssemblyParameters().k = 4;
		List<DirectedEvidence> input = new ArrayList<DirectedEvidence>();
		input.add(SCE(BWD, withSequence("CGTAACCGGTTC", Read(0, 9996, "5S4M2I1M"))));
		ArrayList<SAMRecord> r = Lists.newArrayList(new PositionalAssembler(pc, aes, new SequentialIdGenerator("asm"), input.iterator(), null, null));
		assertEquals(1, r.size());
		assertEquals(10000, r.get(0).getAlignmentEnd());
	}
	@Test
	public void should_assemble_each_read_alignment_only_once() {
		ProcessingContext pc = getContext();
		AssemblyEvidenceSource aes = AES(pc);
		pc.getAssemblyParameters().k = 4;
		String seq = "AACCGGTTAA";
		List<DirectedEvidence> input = new ArrayList<DirectedEvidence>();
		input.add(SCE(FWD, withSequence(seq, Read(0, 10, "5M5S"))[0]));
		input.add(SCE(FWD, withSequence(seq, Read(0, 10, "5M5S"))[0]));
		input.add(SCE(FWD, withSequence(seq, Read(0, 10, "5M5S"))[0]));
		input.add(IE(withSequence(seq, Read(0, 10, "5M100D5M"))[0]));
		input.add(IE(withSequence(seq, Read(0, 10, "5M100D5M"))[0]));
		input.add(IE(withSequence(seq, Read(0, 10, "5M100D5M"))[0]));
		input.sort(DirectedEvidenceOrder.ByStartEndStart2End2);
		SAMRecord r = Lists.newArrayList(new PositionalAssembler(pc, aes, new SequentialIdGenerator("asm"), input.iterator(), null, null)).get(0);
		assertEquals(2, new AssemblyAttributes(r).getSupportingReadCount(6, null, null, null));
	}
	@Test
	public void should_set_assembly_direction() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		pc.getAssemblyParameters().maxExpectedBreakendLengthMultiple = 10;
		List<DirectedEvidence> e = new ArrayList<>();
		// ACGTACGTACGTACGTACGT
		// 12345678901234567890
		// MMMM----MMMM
		//         MMMM>>>>
		e.add(SCE(FWD, withSequence("ACGTCCGGACGG", Read(1, 1, "4M8S"))[0]));
		e.add(SCE(FWD, withSequence("ACGGTTGT", Read(1, 9, "4M4S"))[0]));
		e.sort(DirectedEvidenceOrder.ByStartEndStart2End2);
		ArrayList<SAMRecord> output = Lists.newArrayList(new PositionalAssembler(pc, AES(pc), new SequentialIdGenerator("asm"), e.iterator(), BreakendDirection.Forward, null, null));
		for (SAMRecord r : output) {
			assertEquals(BreakendDirection.Forward, new AssemblyAttributes(r).getAssemblyDirection());
		}
	}
	@Test
	public void should_error_correct_and_downsample_at_complexity_threshold() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().errorCorrection.kmerErrorCorrectionMultiple = 2;
		pc.getAssemblyParameters().positional.maximumNodeDensity = 0.05f;
		String seq = S(RANDOM).substring(1, 51) + S(RANDOM).substring(100, 150);
		List<String> sequences = allSequencesWithinEditDistance(seq, 2);
		AtomicInteger index = new AtomicInteger();
		MockSAMEvidenceSource ses = SES();
		List<DirectedEvidence> e = Lists.newArrayList(sequences.stream()
				.map(s -> (DirectedEvidence)SCE(FWD, ses, withName(s + index.incrementAndGet(), withSequence(s, Read(2, 1, "50M50S")))[0]))
				.collect(Collectors.toList()));
		for (int i = 1; i < 200; i++) {
			e.add(SCE(FWD, ses, Read(2, i * 100, "50M50S")));
		}
		IntervalBed excluded = new IntervalBed(pc.getLinear());
		AssemblyEvidenceSource aes = AES(pc);
		ArrayList<SAMRecord> output = Lists.newArrayList(new PositionalAssembler(pc, aes, new SequentialIdGenerator("asm"), e.iterator(), BreakendDirection.Forward, excluded, null));
		Assert.assertNotEquals(0, excluded.size()); // we should have hit the threshold
		List<Collection<String>> assembledReads = output.stream()
				.map(a -> new AssemblyAttributes(a).getEvidenceIDs(null, null, null, aes))
				.distinct()
				.collect(Collectors.toList());
		Assert.assertNotEquals(sequences.stream().distinct().count(), assembledReads.size()); // make sure we assembled fewer reads
	}
}
