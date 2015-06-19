package au.edu.wehi.idsv.debruijn.positional;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import static org.junit.Assert.*;

import org.junit.Test;

import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.TestHelper;

import com.google.common.collect.Lists;


public class PathCollapseIteratorTest extends TestHelper {
	public static Iterator<KmerPathNode> toPathNodeIterator(List<DirectedEvidence> input, int k, int maxWidth, int maxLength) {
		return new PathNodeIterator(new AggregateNodeIterator(new SupportNodeIterator(k, input.iterator(), maxLength)), maxWidth, maxLength, k);
	}
	public static List<KmerPathNode> toPathNode(List<DirectedEvidence> input, int k, int maxWidth, int maxLength) {
		List<KmerNode> snList = Lists.newArrayList(new SupportNodeIterator(k, input.iterator(), maxLength));
		List<KmerNode> anList = Lists.newArrayList(new AggregateNodeIterator(snList.iterator()));
		List<KmerPathNode> pnList = Lists.newArrayList( new PathNodeIterator(anList.iterator(), maxWidth, maxLength, k));
		assertSameNodes(anList, pnList);
		return pnList;
	}
	private static List<KmerPathNode> go(List<DirectedEvidence> input, int k, int maxWidth, int maxLength, int maxPathCollapseLength, int maxBasesMismatch) {
		List<KmerPathNode> pnList = toPathNode(input, k, maxWidth, maxLength);
		PathCollapseIterator pcit = new PathCollapseIterator(pnList.iterator(), k, maxWidth, maxLength, maxPathCollapseLength, maxBasesMismatch);
		ArrayList<KmerPathNode> result = Lists.newArrayList(pcit);
		assertEquals(totalWeight(pnList), totalWeight(result));
		return result;
	}
	@Test
	public void should_collapse_leaf() {
		int k = 25;
		List<DirectedEvidence> input = new ArrayList<DirectedEvidence>();
		input.add(SCE(FWD, withSequence(S(RANDOM).substring(0, 75), Read(0, 1, "50M25S"))));
		input.add(SCE(FWD, withSequence(S(RANDOM).substring(0, 70) + "AA", Read(0, 1, "50M22S"))));
		
		List<KmerPathNode> result = go(input, k, 1, 200, 5, 2);
		assertEquals(1, result.size());
	}
}
