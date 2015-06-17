package au.edu.wehi.idsv.debruijn.positional;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.junit.Test;

import com.google.common.collect.Lists;

import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.SAMEvidenceSource;
import au.edu.wehi.idsv.TestHelper;


public class SupportNodeIteratorTest extends TestHelper {
	@Test
	public void should_return_kmernodes_in_start_position_order() {
		int k = 4;
		SAMEvidenceSource ses = SES(30, 60);
		List<DirectedEvidence> input = new ArrayList<DirectedEvidence>();
		for (int i = 1; i <= 100; i++) {
			input.add(SCE(FWD, ses, withQual(new byte[] { 0,1,2,3,4,5,6,7,8,9,10}, withSequence("ACGTTATACCG", Read(0, i, "1S4M6S")))));
			input.add(SCE(BWD, ses, withQual(new byte[] { 0,1,2,3,4,5,6,7,8,9,10}, withSequence("ACGTTATACCG", Read(0, i, "1S4M6S")))));
			input.add(NRRP(ses, withQual(new byte[] { 0,1,2,3,4,5,6,7,8,9,10}, withSequence("ACGTTATACCG", DP(0, 21, "1S9M1S", false, 1, 1, "11M", true)))));
			input.add(NRRP(ses, withQual(new byte[] { 0,1,2,3,4,5,6,7,8,9,10}, withSequence("ACGTTATACCG", DP(0, 21, "1S9M1S", false, 1, 1, "11M", false)))));
		}
		Collections.sort(input, DirectedEvidence.ByStartEnd);
		List<KmerSupportNode> output = Lists.newArrayList(new SupportNodeIterator(k, input.iterator(), ses.getMaxConcordantFragmentSize()));
		assertTrue(KmerNode.ByStartPosition.isOrdered(output));
		assertEquals(100 * ( 10-3 + 5-3 + 11-3 + 11-3 ), output.size());
	}
}
