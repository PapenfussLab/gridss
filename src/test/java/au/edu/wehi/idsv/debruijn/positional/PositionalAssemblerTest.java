package au.edu.wehi.idsv.debruijn.positional;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;

import au.edu.wehi.idsv.sam.SamTags;
import com.google.common.collect.Range;
import org.junit.Test;

import com.google.common.collect.Lists;

import au.edu.wehi.idsv.AssemblyAttributes;
import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.BreakendSummary;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.DirectedEvidenceOrder;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.SequentialIdGenerator;
import au.edu.wehi.idsv.SingleReadEvidence;
import au.edu.wehi.idsv.TestHelper;
import htsjdk.samtools.SAMRecord;


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
		input.sort(DirectedEvidenceOrder.ByStartEnd);
		ArrayList<SAMRecord> r = Lists.newArrayList(new PositionalAssembler(pc, aes, new SequentialIdGenerator("asm"), input.iterator()));
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
		input.sort(DirectedEvidenceOrder.ByStartEnd);
		List<SingleReadEvidence> r = asAssemblyEvidence(aes, Lists.newArrayList(new PositionalAssembler(pc, aes, new SequentialIdGenerator("asm"), input.iterator())));
		assertEquals(1, r.size());
		assertEquals(new BreakendSummary(0, FWD, 14), r.get(0).getBreakendSummary());
		// anchor length to match breakend length
		assertEquals("AACGTT", S(r.get(0).getAnchorSequence()));
		assertEquals("GGTTAA", S(r.get(0).getBreakendSequence()));
		assertEquals("AACGTTGGTTAA", S(r.get(0).getSAMRecord().getReadBases()));
		int[] scCount = r.get(0).getSAMRecord().getSignedIntArrayAttribute(SamTags.ASSEMBLY_SOFTCLIP_COUNT);
		assertArrayEquals(new int[] {
			  // A A C G T T G G T T A A
				1,1,2,2,2,2,2,2,2,2,2,1,0, // category 0
				0,0,0,0,0,0,0,0,0,0,0,0,0, // category 1
			}, scCount);
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
		input.sort(DirectedEvidenceOrder.ByStartEnd);
		List<SingleReadEvidence> r = asAssemblyEvidence(aes, Lists.newArrayList(new PositionalAssembler(pc, aes, new SequentialIdGenerator("asm"), input.iterator())));
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
		ArrayList<SAMRecord> r = Lists.newArrayList(new PositionalAssembler(pc, aes, new SequentialIdGenerator("asm"), input.iterator()));
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
		ArrayList<SAMRecord> r = Lists.newArrayList(new PositionalAssembler(pc, aes, new SequentialIdGenerator("asm"), input.iterator()));
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
		input.sort(DirectedEvidenceOrder.ByStartEnd);
		SAMRecord r = Lists.newArrayList(new PositionalAssembler(pc, aes, new SequentialIdGenerator("asm"), input.iterator())).get(0);
		assertEquals(2, (int)new AssemblyAttributes(r).getSupportingReadCount(Range.closedOpen(0, r.getReadLength()), null, null, Math::max).getRight());
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
		e.sort(DirectedEvidenceOrder.ByStartEnd);
		ArrayList<SAMRecord> output = Lists.newArrayList(new PositionalAssembler(pc, AES(pc), new SequentialIdGenerator("asm"), e.iterator(), BreakendDirection.Forward));
		for (SAMRecord r : output) {
			assertEquals(BreakendDirection.Forward, new AssemblyAttributes(r).getAssemblyDirection());
		}
	}
}
