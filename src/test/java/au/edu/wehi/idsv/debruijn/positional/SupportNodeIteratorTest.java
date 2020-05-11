package au.edu.wehi.idsv.debruijn.positional;

import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.SAMEvidenceSource;
import au.edu.wehi.idsv.TestHelper;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;


public class SupportNodeIteratorTest extends TestHelper {
	public static List<DirectedEvidence> scrp(int k, String sequence, int minFragSize, int maxFragSize) {
		SAMEvidenceSource ses = SES(minFragSize, maxFragSize);
		byte[] qual = new byte[sequence.length()];
		for (int i = 0; i < qual.length; i++) {
			qual[i] = (byte)i;
		}
		List<DirectedEvidence> input = new ArrayList<DirectedEvidence>();
		for (int i = 1; i <= 100; i++) {
			input.add(SCE(FWD, ses, withQual(qual, withSequence(sequence, Read(0, i, "1S4M6S")))));
			input.add(SCE(BWD, ses, withQual(qual, withSequence(sequence, Read(0, i, "1S4M6S")))));
			input.add(NRRP(ses, withQual(qual, withSequence(sequence, DP(0, i, "1S9M1S", false, 1, 1, "11M", true)))));
			input.add(NRRP(ses, withQual(qual, withSequence(sequence, DP(0, i, "1S9M1S", false, 1, 1, "11M", false)))));
		}
		Collections.sort(input, DirectedEvidence.ByStartEnd);
		return input;
	}
	@Test
	public void should_return_kmernodes_in_start_position_order() {
		int k = 4;
		List<DirectedEvidence> input = scrp(k, "ACGTTATACCG", 30, 60);
		List<KmerSupportNode> output = Lists.newArrayList(new SupportNodeIterator(k, input.iterator(), 60, null, false, 0));
		assertTrue(KmerNodeUtil.ByLastStart.isOrdered(output));
		assertEquals(100 * ( 10-3 + 5-3 + 11-3 + 11-3 ), output.size());
	}
	@Test
	public void includePairAnchors_should_determine_rp_anchor_inclusion() {
		int k = 4;
		List<DirectedEvidence> input = ImmutableList.of(NRRP(OEA(0, 1, "10M", true)));
		assertEquals(1 * (10-3), Lists.newArrayList(new SupportNodeIterator(k, input.iterator(), 60, null, false, 0)).size());
		assertEquals(2 * (10-3), Lists.newArrayList(new SupportNodeIterator(k, input.iterator(), 60, null, true, 0)).size());
	}
	@Test
	public void should_only_include_ungapped_anchor_subalignment_kmers() {
		int k = 4;
		List<DirectedEvidence> input;
		input = ImmutableList.of(NRRP(OEA(0, 1, "6M4S", true)));
		assertEquals((10-3) + (6-3), Lists.newArrayList(new SupportNodeIterator(k, input.iterator(), 60, null, true, 0)).size());
		input = ImmutableList.of(NRRP(OEA(0, 1, "4M2I4M", true)));
		assertEquals((10-3) + 1, Lists.newArrayList(new SupportNodeIterator(k, input.iterator(), 60, null, true, 0)).size());
		input = ImmutableList.of(NRRP(OEA(0, 1, "5M1D5M", true)));
		assertEquals((10-3) + 2, Lists.newArrayList(new SupportNodeIterator(k, input.iterator(), 60, null, true, 0)).size());
	}
	@Test
	public void should_not_include_anchor_base_mismatches_with_5bp() {
		int k = 4;
		List<DirectedEvidence> input;
		input = ImmutableList.of(NRRP(withSequence("AAAAAAAAAT", OEA(0, 1, "10M", true))));
		assertEquals((10-3) + (10-3-1), Lists.newArrayList(new SupportNodeIterator(k, input.iterator(), 60, null, true, 5)).size());
		input = ImmutableList.of(NRRP(withSequence("AAAAAAAATT", OEA(0, 1, "10M", true))));
		assertEquals((10-3) + (10-3-2), Lists.newArrayList(new SupportNodeIterator(k, input.iterator(), 60, null, true, 5)).size());
		input = ImmutableList.of(NRRP(withSequence("AAAAAAAATA", OEA(0, 1, "10M", true))));
		assertEquals((10-3) + (10-3-2), Lists.newArrayList(new SupportNodeIterator(k, input.iterator(), 60, null, true, 5)).size());
	}
}
