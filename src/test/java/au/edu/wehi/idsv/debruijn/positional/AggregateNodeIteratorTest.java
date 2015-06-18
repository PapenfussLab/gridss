package au.edu.wehi.idsv.debruijn.positional;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.junit.Test;

import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.SAMEvidenceSource;
import au.edu.wehi.idsv.TestHelper;

import com.google.common.collect.Lists;


public class AggregateNodeIteratorTest extends TestHelper {
	@Test
	public void should_aggregate_kmer_weights() {
		int k = 4;
		SAMEvidenceSource ses = SES(30, 60);
		List<DirectedEvidence> input = new ArrayList<DirectedEvidence>();
		input.add(NRRP(ses, withQual(new byte[] { 1,1,1,1,1}, withSequence("ACGTT", DP(0, 1, "5M", false, 1, 1, "5M", false)))));
		input.add(NRRP(ses, withQual(new byte[] { 1,1,1,1,1}, withSequence("ACGTT", DP(0, 2, "5M", false, 1, 1, "5M", false)))));
		Collections.sort(input, DirectedEvidence.ByStartEnd);
		List<KmerSupportNode> snList = Lists.newArrayList(new SupportNodeIterator(k, input.iterator(), ses.getMaxConcordantFragmentSize()));
		List<KmerAggregateNode> anList = Lists.newArrayList(new AggregateNodeIterator(snList.iterator()));
		assertTrue(KmerNode.ByStartPosition.isOrdered(anList));
		assertEquals(4, anList.size());
		assertEquals(6, anList.size());
	}
	@Test
	public void total_weight_should_remain_constant() {
		List<KmerSupportNode> snList = Lists.newArrayList(new SupportNodeIterator(4, SupportNodeIteratorTest.scrp(4, "ACGTTATACCG", 30, 60).iterator(), 60));
		List<KmerAggregateNode> anList = Lists.newArrayList(new AggregateNodeIterator(snList.iterator()));
		assertEquals(
				snList.stream().mapToInt(n -> (n.endPosition() - n.startPosition() + 1) * n.weight()).sum(),
				anList.stream().mapToInt(n -> (n.endPosition() - n.startPosition() + 1) * n.weight()).sum());
	}
}
