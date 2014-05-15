package au.edu.wehi.socrates.debruijn;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import java.util.List;

import org.junit.Test;

import au.edu.wehi.socrates.BreakendDirection;
import au.edu.wehi.socrates.DirectedBreakpointAssembly;
import au.edu.wehi.socrates.DirectedEvidence;
import au.edu.wehi.socrates.SoftClipEvidence;
import au.edu.wehi.socrates.TestHelper;

import com.google.common.collect.Lists;


public class DeBruijnAssemblerTest extends TestHelper {
	public List<DirectedBreakpointAssembly> go(int k, DirectedEvidence... evidence) {
		List<DirectedBreakpointAssembly> list = Lists.newArrayList();
		DeBruijnAssembler assembler = new DeBruijnAssembler(getContext(), k);
		for (DirectedEvidence e : evidence) {
			Iterable<DirectedBreakpointAssembly> it = assembler.addEvidence(e);
			if (it != null) {
				for (DirectedBreakpointAssembly ass : it) {
					list.add(ass);
				}
			}
		}
		Iterable<DirectedBreakpointAssembly> it = assembler.endOfEvidence();
		if (it != null) {
			for (DirectedBreakpointAssembly ass : it) {
				list.add(ass);
			}
		}
		return list;
	}
	@Test
	public void should_not_call_if_no_evidence() {
		List<DirectedBreakpointAssembly> r = go(3, new SoftClipEvidence(getContext(), BreakendDirection.Backward, Read(1, 1, "2S5M")));
		assertEquals(0, r.size());
	}
	@Test
	public void should_not_call_single_soft_clip() {
		List<DirectedBreakpointAssembly> r = go(3,
				SCE(BreakendDirection.Forward, withQual(new byte[] { 5,5,5,5,5,5 }, withSequence("AACGTG", Read(0, 1, "1M5S"))))
		); 
		assertEquals(0, r.size());
	}
	@Test
	public void should_not_call_if_no_breakpoint_assembly() {
		// polyA reads assemble as anchor kmer
		List<DirectedBreakpointAssembly> r = go(3,
				new SoftClipEvidence(getContext(), BreakendDirection.Backward, Read(0, 1, "2S5M")),
				new SoftClipEvidence(getContext(), BreakendDirection.Backward, Read(0, 1, "3S5M"))
		); 
		assertEquals(0, r.size());
	}
	@Test
	public void should_not_call_unanchored_evidence() {
		// polyA reads assemble as anchor kmer
		List<DirectedBreakpointAssembly> r = go(3,
				NRRP(withSequence("CATG", OEA(1, 1, "10M", true))),
				NRRP(withSequence("CATGAT", OEA(1, 1, "10M", true)))
		);
		assertEquals(0, r.size());
	}
	@Test
	public void should_call_multiple_soft_clips() {
		// polyA reads assemble as anchor kmer
		List<DirectedBreakpointAssembly> r = go(3,
				SCE(BreakendDirection.Forward, withQual(new byte[] { 5,5,5,5,5,5 }, withSequence("AACGTG", Read(0, 1, "1M5S")))),
				SCE(BreakendDirection.Forward, withQual(new byte[] { 5,5,5,5,5,5,5 }, withSequence("AACGTGA", Read(0, 1, "1M6S"))))
		); 
		assertEquals(1, r.size());
		assertEquals("ACGTGA", r.get(0).getBreakpointSequenceString());
		assertEquals("A", r.get(0).getAnchorSequenceString());
		assertEquals(2, r.get(0).getConsensusReadCount());
		assertEquals(BreakendDirection.Forward, r.get(0).getBreakendSummary().direction);
		assertEquals(1, r.get(0).getBreakendSummary().start);
		assertEquals(1, r.get(0).getBreakendSummary().end);
		assertEquals(0, r.get(0).getBreakendSummary().referenceIndex);
		assertEquals(7, r.get(0).getBreakendSummary().qual, 0);
	}
	//@Test
	public void should_call_with_breakpoint_quality() {
		List<DirectedBreakpointAssembly> r = go(3,
				SCE(BreakendDirection.Forward, withQual(new byte[] { 5,5,5,5,5,5 }, withSequence("AACGTG", Read(0, 1, "1M5S")))),
				SCE(BreakendDirection.Forward, withQual(new byte[] { 5,5,5,5,5,5,5 }, withSequence("AACGTGA", Read(0, 1, "1M6S"))))
		); 
		assertArrayEquals(new byte[] { 10,10,10,10,10,5 }, r.get(0).getBreakpointQuality());
	}
}
