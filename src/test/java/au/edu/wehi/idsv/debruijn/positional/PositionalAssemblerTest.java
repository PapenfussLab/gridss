package au.edu.wehi.idsv.debruijn.positional;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.BreakendSummary;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.DirectedEvidenceOrder;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.SAMRecordAssemblyEvidence;
import au.edu.wehi.idsv.TestHelper;

import com.google.common.collect.Lists;


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
		ArrayList<SAMRecordAssemblyEvidence> r = Lists.newArrayList(new PositionalAssembler(pc, aes, input.iterator()));
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
		input.add(SCE(FWD, withSequence("ACGTTGGTTA", Read(0, 10, "5M5S"))[0]));
		input.add(SCE(FWD, withSequence("TTTTTGCAACGTTGGTTAA", Read(0, 2, "13M6S"))[0]));
		input.sort(DirectedEvidenceOrder.ByStartEnd);
		ArrayList<SAMRecordAssemblyEvidence> r = Lists.newArrayList(new PositionalAssembler(pc, aes, input.iterator()));
		assertEquals(1, r.size());
		assertEquals(new BreakendSummary(0, FWD, 14, 14), r.get(0).getBreakendSummary());
		// anchor length to match breakend length
		assertEquals("AACGTT", S(r.get(0).getAssemblyAnchorSequence()));
		assertEquals("GGTTAA", S(r.get(0).getBreakendSequence()));
		assertEquals("AACGTTGGTTAA", S(r.get(0).getAssemblySequence()));
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
		ArrayList<SAMRecordAssemblyEvidence> r = Lists.newArrayList(new PositionalAssembler(pc, aes, input.iterator()));
		assertEquals(2, r.size());
		// race condition w.r.t which assembly returns first
		assertEquals(new BreakendSummary(0, FWD, 14, 14), r.get(0).getBreakendSummary());
		// anchor length to match breakend length
		assertEquals("AACGTT", S(r.get(0).getAssemblyAnchorSequence()));
		assertEquals("GGTTAA", S(r.get(0).getBreakendSequence()));
		assertEquals("AACGTTGGTTAA", S(r.get(0).getAssemblySequence()));
	}
}
