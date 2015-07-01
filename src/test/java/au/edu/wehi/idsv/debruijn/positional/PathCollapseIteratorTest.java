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
	public static List<KmerPathNode> simplify(List<KmerPathNode> list) {
		ArrayList<KmerPathNode> result = Lists.newArrayList(new PathSimplificationIterator(list.iterator(), 10000, 10000));
		result.sort(KmerNodeUtil.ByFirstStartEndKmerReference);
		return result;
	}
	public static Iterator<KmerPathNode> toPathNodeIterator(List<DirectedEvidence> input, int k, int maxWidth, int maxLength) {
		return new PathNodeIterator(new AggregateNodeIterator(new SupportNodeIterator(k, input.iterator(), 10000)), maxLength, k);
	}
	public static List<KmerPathNode> toPathNode(List<DirectedEvidence> input, int k, int maxWidth, int maxLength) {
		List<KmerNode> snList = Lists.newArrayList(new SupportNodeIterator(k, input.iterator(), 10000));
		List<KmerNode> anList = Lists.newArrayList(new AggregateNodeIterator(snList.iterator()));
		List<KmerPathNode> pnList = Lists.newArrayList(new PathNodeIterator(anList.iterator(), maxWidth, k));
		assertSameNodes(anList, pnList);
		assert(KmerNodeUtil.ByFirstStartKmer.isOrdered(pnList));
		return pnList;
	}
	private static List<KmerPathNode> go(List<DirectedEvidence> input, int k, int maxLength, int maxPathCollapseLength, int maxBasesMismatch) {
		List<KmerPathNode> pnList = asCheckedKPN(k, maxLength, input.toArray(new DirectedEvidence[0]));
		return go(k, maxPathCollapseLength, maxBasesMismatch, pnList);
	}
	private static List<KmerPathNode> go(int k, int maxPathCollapseLength, int maxBasesMismatch, List<KmerPathNode> pnList) {
		int pnTotalWeight = totalWeight(pnList);
		PathCollapseIterator pcit = new PathCollapseIterator(pnList.iterator(), k, maxPathCollapseLength, maxBasesMismatch, false);
		ArrayList<KmerPathNode> result = Lists.newArrayList(pcit);
		assertEquals(pnTotalWeight, totalWeight(result));
		return result;
	}
	@Test
	public void should_collapse_path() {
		int k = 4;
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(k, "GTAC", 1, 10, false)); 
		input.add(KPN(k, "TACTAAA", 3, 11, false, 2));
		input.add(KPN(k, "TACGAAA", 4, 5, false));
		input.add(KPN(k, "AAAT", 6, 15, false));
		KmerPathNode.addEdge(input.get(0), input.get(1));
		KmerPathNode.addEdge(input.get(0), input.get(2));
		KmerPathNode.addEdge(input.get(1), input.get(3));
		KmerPathNode.addEdge(input.get(2), input.get(3));
		List<KmerPathNode> result = simplify(go(k, 100, 100, input));
		assertTrue(result.contains(KPN(k, "GTAC", 1, 10, false)));
		assertTrue(result.contains(KPN(k, "TACTAAA", 3, 3, false, 2)));
		assertTrue(result.contains(KPN(k, "TACTAAA", 4, 5, false, 3)));
		assertTrue(result.contains(KPN(k, "TACTAAA", 6, 11, false, 2)));
		assertTrue(result.contains(KPN(k, "AAAT", 6, 15, false)));
	}
	@Test
	public void collapse_should_not_occur_above_maxBasesMismatch() {
		int k = 25;
		List<DirectedEvidence> input = new ArrayList<DirectedEvidence>();
		input.add(SCE(FWD, withSequence(S(RANDOM).substring(0, 75), Read(0, 1, "50M25S"))));
		input.add(SCE(FWD, withSequence(S(RANDOM).substring(0, 70) + "GG", Read(0, 1, "50M22S"))));
		
		List<KmerPathNode> result = go(input, k, 200, 200, 1);
		assertEquals(4, result.size());
	}
	@Test
	public void collapse_should_not_occur_above_min_path_length() {
		int k = 25;
		List<DirectedEvidence> input = new ArrayList<DirectedEvidence>();
		input.add(SCE(FWD, withSequence(S(RANDOM).substring(0, 75), Read(0, 1, "50M25S"))));
		input.add(SCE(FWD, withSequence(S(RANDOM).substring(0, 70) + "GG", Read(0, 1, "50M22S"))));
		
		List<KmerPathNode> result = go(input, k, 200, 1, 100);
		assertEquals(4, result.size());
	}
	@Test
	public void should_split_position_mid() {
		int k = 4;
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(k, "GTACC", 1, 10, false));
		input.add(KPN(k, "ACCTT", 2, 11, false)); // this gets split on position as weight increases at position 5
		input.add(KPN(k, "ACCG", 5, 5, false));
		KmerPathNode.addEdge(input.get(0), input.get(1));
		KmerPathNode.addEdge(input.get(0), input.get(2));
		List<KmerPathNode> result = go(k, 5, 1, input);
		assertSame(result.get(0), KPN(k, "GTACC", 1, 10, false));
		assertSame(result.get(1), KPN(k, "ACCTT", 2, 4, false));
		assertSame(result.get(2), KPN(k, "ACCT", 5, 5, false, new int[] { 2, }));
		assertSame(result.get(3), KPN(k, "CCTT", 6, 6, false));
		assertSame(result.get(4), KPN(k, "ACCTT", 6, 11, false));
	}
	@Test
	public void should_split_position_start() {
		int k = 4;
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(k, "GTACC", 0, 9, false));
		input.add(KPN(k, "ACCTT", 2, 11, false)); // this gets split on position as weight increases at position 5
		input.add(KPN(k, "ACCG", 2, 5, false));
		KmerPathNode.addEdge(input.get(0), input.get(1));
		KmerPathNode.addEdge(input.get(0), input.get(2));
		List<KmerPathNode> result = go(k, 5, 1, input);
		assertEquals(4, result.size());
		assertSame(result.get(0), KPN(k, "GTACC", 0, 9, false));
		assertSame(result.get(1), KPN(k, "ACCT", 2, 5, false, new int[] { 2 }));
		assertSame(result.get(2), KPN(k, "CCTT", 3, 6, false));
		assertSame(result.get(3), KPN(k, "ACCTT", 6, 11, false));
	}
	@Test
	public void should_split_position_end() {
		int k = 4;
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(k, "GTACC", 1, 10, false));
		input.add(KPN(k, "ACCT", 2, 11, false)); // this gets split on position as weight increases at position 5
		input.add(KPN(k, "ACCGG", 5, 11, false));
		KmerPathNode.addEdge(input.get(0), input.get(1));
		KmerPathNode.addEdge(input.get(0), input.get(2));
		List<KmerPathNode> result = go(k, 10, 1, input);
		assertEquals(4, result.size());
		assertSame(result.get(0), KPN(k, "GTACC", 1, 10, false));
		assertSame(result.get(1), KPN(k, "ACCT", 2, 4, false));
		assertSame(result.get(2), KPN(k, "ACCG", 5, 11, false, new int[] { 2}));
		assertSame(result.get(3), KPN(k, "CCGG", 6, 12, false));
	}
	@Test
	public void should_traverse_forks() {
		//                  G - C
		//                 /  \
		// AAAA - TTGGCCA      T
		//                \
		//                  TAA
		int k = 4;          
		List<DirectedEvidence> input = new ArrayList<DirectedEvidence>();
		input.add(SCE(FWD, withSequence("AAAATTGGCCAGC", Read(0, 1, "4M9S"))));
		input.add(SCE(FWD, withSequence("AAAATTGGCCAGT", Read(0, 1, "4M9S"))));
		input.add(SCE(FWD, withSequence("AAAATTGGCCATAA", Read(0, 1, "4M10S"))));
		List<KmerPathNode> result = go(input, k, 200, 200, 1);
		result = Lists.newArrayList(new PathSimplificationIterator(result.iterator(), 100, 100));
		assertEquals(4, result.size());
	}
	@Test
	public void should_not_collapse_higher_weighted_leaf() {
		int k = 4;   
		List<DirectedEvidence> input = new ArrayList<DirectedEvidence>();
		input.add(SCE(FWD, withSequence("AAAATTGGCCAGCTCATTTT", Read(0, 1, "4M16S"))));
		input.add(SCE(FWD, withSequence("AAAATTGGCCAGCTCAGGGG", Read(0, 1, "4M16S"))));
		input.add(SCE(FWD, withSequence("AAAATTGGCCAGCTCT", Read(0, 1, "4M12S"))));
		input.add(SCE(FWD, withSequence("AAAATTGGCCAGCTCT", Read(0, 1, "4M12S"))));
		input.add(SCE(FWD, withSequence("AAAATTGGCCAGCTCT", Read(0, 1, "4M12S"))));
		input.add(SCE(FWD, withSequence("AAAATTGGCCAGCTCT", Read(0, 1, "4M12S"))));
		input.add(SCE(FWD, withSequence("AAAATTGGCCAGCTCT", Read(0, 1, "4M12S"))));
		input.add(SCE(FWD, withSequence("AAAATTGGCCAGCTCT", Read(0, 1, "4M12S"))));
		input.add(SCE(FWD, withSequence("AAAATTGGCCAGCTCT", Read(0, 1, "4M12S"))));
		input.add(SCE(FWD, withSequence("AAAATTGGCCAGCTCT", Read(0, 1, "4M12S"))));
		input.add(SCE(FWD, withSequence("AAAATTGGCCAGCTCT", Read(0, 1, "4M12S"))));
		List<KmerPathNode> result = go(input, k, 200, 200, 1);
		assertEquals(6, result.size());
	}
	@Test
	public void should_collapse_lower_higher_weighted_leaf() {
		int k = 4;   
		List<DirectedEvidence> input = new ArrayList<DirectedEvidence>();
		input.add(SCE(FWD, withSequence("AAAATTGGCCAGCTCATTTT", Read(0, 1, "4M16S"))));
		input.add(SCE(FWD, withSequence("AAAATTGGCCAGCTCAGGGG", Read(0, 1, "4M16S"))));
		input.add(SCE(FWD, withSequence("AAAATTGGCCAGCTCT", Read(0, 1, "4M12S"))));
		List<KmerPathNode> result = go(input, k, 200, 200, 1);
		assertEquals(5, result.size());
	}
	@Test
	public void should_not_compress_loops() {
		int k = 4;
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(k, "AAAA", 0, 10, false));
		input.add(KPN(k, "AAAG", 0, 10, false));
		input.add(KPN(k, "AAGC", 0, 15, false));
		KmerPathNode.addEdge(input.get(0), input.get(0));
		KmerPathNode.addEdge(input.get(0), input.get(1));
		KmerPathNode.addEdge(input.get(1), input.get(2));
		List<KmerPathNode> result = go(k, 100, 100, input);
		// AAAA->AAAA
		// can collapse with
		// AAAA->AAAG
		// TODO: actual test case that tests loop compression
		assertTrue(result.stream().allMatch(pn -> pn.collapsedKmers().size() == 0));
	}
}
