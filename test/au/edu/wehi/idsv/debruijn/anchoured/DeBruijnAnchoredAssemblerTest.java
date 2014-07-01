package au.edu.wehi.idsv.debruijn.simple;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import java.util.List;

import org.junit.Test;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.SoftClipEvidence;
import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;
import au.edu.wehi.idsv.debruijn.anchoured.DeBruijnAnchoredAssembler;
import au.edu.wehi.idsv.vcf.VcfAttributes;
import au.edu.wehi.idsv.TestHelper;

import com.google.common.collect.Lists;


public class DeBruijnAnchoredAssemblerTest extends TestHelper {
	public List<VariantContextDirectedBreakpoint> go(int k, DirectedEvidence... evidence) {
		List<VariantContextDirectedBreakpoint> list = Lists.newArrayList();
		DeBruijnWindowedAssembler assembler = new DeBruijnWindowedAssembler(getContext(), k);
		for (DirectedEvidence e : evidence) {
			Iterable<VariantContextDirectedBreakpoint> it = assembler.addEvidence(e);
			if (it != null) {
				for (VariantContextDirectedBreakpoint ass : it) {
					list.add(ass);
				}
			}
		}
		Iterable<VariantContextDirectedBreakpoint> it = assembler.endOfEvidence();
		if (it != null) {
			for (VariantContextDirectedBreakpoint ass : it) {
				list.add(ass);
			}
		}
		return list;
	}
	@Test
	public void should_not_call_if_no_evidence() {
		List<VariantContextDirectedBreakpoint> r = go(3, new SoftClipEvidence(getContext(), BreakendDirection.Backward, Read(1, 1, "2S5M")));
		assertEquals(0, r.size());
	}
	@Test
	public void should_not_call_single_soft_clip() {
		List<VariantContextDirectedBreakpoint> r = go(3,
				SCE(BreakendDirection.Forward, withQual(new byte[] { 5,5,5,5,5,5 }, withSequence("AACGTG", Read(0, 1, "1M5S"))))
		); 
		assertEquals(0, r.size());
	}
	@Test
	public void should_not_call_if_no_breakpoint_assembly() {
		// polyA reads assemble as anchor kmer
		List<VariantContextDirectedBreakpoint> r = go(3,
				new SoftClipEvidence(getContext(), BreakendDirection.Backward, Read(0, 1, "2S5M")),
				new SoftClipEvidence(getContext(), BreakendDirection.Backward, Read(0, 1, "3S5M"))
		); 
		assertEquals(0, r.size());
	}
	@Test
	public void should_not_call_unanchored_evidence() {
		// polyA reads assemble as anchor kmer
		List<VariantContextDirectedBreakpoint> r = go(3,
				NRRP(withSequence("CATG", OEA(1, 1, "10M", true))),
				NRRP(withSequence("CATGAT", OEA(1, 1, "10M", true)))
		);
		assertEquals(0, r.size());
	}
	@Test
	public void should_call_multiple_soft_clips() {
		// polyA reads assemble as anchor kmer
		List<VariantContextDirectedBreakpoint> r = go(3,
				SCE(BreakendDirection.Forward, withQual(new byte[] { 5,5,5,5,5,5 }, withSequence("AACGTG", Read(0, 1, "1M5S")))),
				SCE(BreakendDirection.Forward, withQual(new byte[] { 5,5,5,5,5,5,5 }, withSequence("AACGTGA", Read(0, 1, "1M6S"))))
		); 
		assertEquals(1, r.size());
		assertEquals("ACGTGA", r.get(0).getBreakpointSequenceString());
		assertEquals("A", r.get(0).getAnchorSequenceString());
		assertEquals(2, r.get(0).getBreakendSummary().evidence.get(VcfAttributes.ASSEMBLY_READS));
		assertEquals(BreakendDirection.Forward, r.get(0).getBreakendSummary().direction);
		assertEquals(1, r.get(0).getBreakendSummary().start);
		assertEquals(1, r.get(0).getBreakendSummary().end);
		assertEquals(0, r.get(0).getBreakendSummary().referenceIndex);
	}
	@Test
	public void should_call_with_breakpoint_quality() {
		List<VariantContextDirectedBreakpoint> r = go(3,
				SCE(BreakendDirection.Forward, withQual(new byte[] { 1,2,3,4,5,6 }, withSequence("AACGTG", Read(0, 2, "1M5S")))),
				SCE(BreakendDirection.Forward, withQual(new byte[] { 6,7,8,9,10,11,12,13,14 }, withSequence("TAACGTGAT", Read(0, 1, "2M6S"))))
		);
		// kmer qual = sum of min base quals
		// end is padded
		// first two bases are ignored since they're part of the anchor assembly
		assertArrayEquals(new byte[] { /*6,1+7,*/2+8,3+9,4+10,11,12,12,12 }, r.get(0).getBreakpointQuality());
		// TODO: pad both ends so qual is balanced
		// TODO: better base qual weighting
	}
	@Test
	public void id_should_contain_assembler_name_position_direction() {
		List<VariantContextDirectedBreakpoint> r = go(3,
				SCE(BreakendDirection.Forward, withQual(new byte[] { 5,5,5,5,5,5 }, withSequence("AACGTG", Read(0, 1, "1M5S")))),
				SCE(BreakendDirection.Forward, withQual(new byte[] { 5,5,5,5,5,5,5 }, withSequence("AACGTGA", Read(0, 1, "1M6S"))))
		); 
		assertEquals("debruijn-polyA:1-f", r.get(0).getEvidenceID());
	}
}
