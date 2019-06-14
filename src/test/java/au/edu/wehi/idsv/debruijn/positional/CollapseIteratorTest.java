package au.edu.wehi.idsv.debruijn.positional;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.junit.Test;

import com.google.common.collect.Lists;

import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.TestHelper;


public abstract class CollapseIteratorTest extends TestHelper {
	public static List<KmerPathNode> simplify(List<KmerPathNode> list) {
		ArrayList<KmerPathNode> result = Lists.newArrayList(new PathSimplificationIterator(list.iterator(), 10000, 10000));
		result.sort(KmerNodeUtil.ByFirstStartEndKmerReference);
		return result;
	}
	public static Iterator<KmerPathNode> toPathNodeIterator(List<DirectedEvidence> input, int k, int maxWidth, int maxLength) {
		return new PathNodeIterator(new AggregateNodeIterator(new SupportNodeIterator(k, input.iterator(), 10000, null, false, 0)), maxLength, k);
	}
	public static List<KmerPathNode> toPathNode(List<DirectedEvidence> input, int k, int maxWidth, int maxLength) {
		List<KmerNode> snList = Lists.newArrayList(new SupportNodeIterator(k, input.iterator(), 10000, null, false, 0));
		List<KmerNode> anList = Lists.newArrayList(new AggregateNodeIterator(snList.iterator()));
		List<KmerPathNode> pnList = Lists.newArrayList(new PathNodeIterator(anList.iterator(), maxWidth, k));
		assertSameNodes(anList, pnList);
		assert(KmerNodeUtil.ByFirstStartKmer.isOrdered(pnList));
		return pnList;
	}
	protected List<KmerPathNode> go(List<DirectedEvidence> input, int k, int maxLength, int maxPathCollapseLength, int maxBasesMismatch) {
		List<KmerPathNode> pnList = asCheckedKPN(k, maxLength, input.toArray(new DirectedEvidence[0]));
		return go(k, maxPathCollapseLength, maxBasesMismatch, pnList);
	}
	protected abstract CollapseIterator create(
			Iterator<KmerPathNode> it,
			int k,
			int maxPathCollapseLength,
			int maxBasesMismatch,
			boolean bubblesAndLeavesOnly,
			double minimumPathNodeEntropy);
	protected List<KmerPathNode> go(int k, int maxPathCollapseLength, int maxBasesMismatch, List<KmerPathNode> pnList) {
		assertDisjoint(pnList); // precondition: is the test case well formed
		int pnTotalWeight = totalWeight(pnList);
		CollapseIterator pcit = create(pnList.iterator(), k, maxPathCollapseLength, maxBasesMismatch, false, 0);
		ArrayList<KmerPathNode> result = Lists.newArrayList(pcit);
		assertEquals(pnTotalWeight, totalWeight(result));
		for (int i = 1; i < result.size(); i++) {
			assertTrue(result.get(i - 1).firstStart() <= result.get(i).firstStart());
		}
		result.sort(KmerNodeUtil.ByFirstStartEndKmerReference);
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
	@Test
	public void should_collapse_kmers() {
		// mostly this test is about checking assertions
		for (int k = 1; k <= 4; k++) {
			List<DirectedEvidence> input = new ArrayList<DirectedEvidence>();
			for (long encoded = 0; encoded < 1 << (2 * k); encoded++) {
				String seq = K(k, encoded);
				input.add(NRRP(withSequence(seq, DP(0, 1, String.format("%dM", k), false, 1, 1, String.format("%dM", k), true))));
			}
			go(input, k, 100, 100, 4);
		}
	}
	@Test
	public void should_collapse_bubbles() {
		int k = 10;
		int sclen = 20;
		int anchorlen = 100;
		boolean twoBasesDifferent = false;
		String seq = S(RANDOM).substring(0, sclen + anchorlen);
		List<DirectedEvidence> input = new ArrayList<DirectedEvidence>();
		for (int i = 0; i < sclen; i++) {
			for (int j = 0; j < sclen; j++) {
				for (byte ci : new byte[] { 'a', 'c', 'g', 't' }) {
					for (byte cj : new byte[] { 'a', 'c', 'g', 't' }) {
						byte[] readSeq = B(seq);
						readSeq[100 + i] = ci;
						if (twoBasesDifferent) readSeq[100 + j] = cj;
						input.add(SCE(FWD, withSequence(readSeq, Read(0, 1, String.format("%dM%dS", sclen, anchorlen)))));
					}
				}
			}
		}
		// should all have collapsed since all initial paths are 2 bases from the start
		// unfortunately, the order of path collapse matters and path collapse happens
		// as soon as a pair of eligible paths are found.
		// To fix this would be computationally expensive as all paths would need
		// to be enumerated before collapse
		go(input, k, 200, 200, 4);
		// assertEquals(2, result.size());
	}
	@Test
	public void should_collapse_simple_leaf_fwd() {
		int k = 4;
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(k, "AAATT", 1, 10, false));
		input.add(KPN(k, "ATTGC", 3, 12, false));
		input.add(KPN(k, "ATTCC", 3, 12, false));
		KmerPathNode.addEdge(input.get(0), input.get(1));
		KmerPathNode.addEdge(input.get(0), input.get(2));
		List<KmerPathNode> result = simplify(go(k, 100, 100, input));
		assertEquals(1, result.size());
	}
	@Test
	public void should_collapse_simple_leaf_bwd() {
		int k = 4;
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(k, "AAATT", 1, 10, false));
		input.add(KPN(k, "ATTGC", 3, 12, false));
		input.add(KPN(k, "CGATT", 1, 10, false));
		KmerPathNode.addEdge(input.get(0), input.get(1));
		KmerPathNode.addEdge(input.get(2), input.get(1));
		List<KmerPathNode> result = simplify(go(k, 100, 100, input));
		assertEquals(1, result.size());
	}
	@Test
	public void should_collapse_multi_path_leaf_fwd() {
		int k = 4;
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(k, "AGATT", 1, 3, false));
		input.add(KPN(k, "ATTCCACGTACGT", 3, 5, false));
		input.add(KPN(k, "ATTGG", 3, 5, false));
		input.add(KPN(k, "TGGAA", 5, 7, false));
		input.add(KPN(k, "GAATT", 7, 9, false));
		KmerPathNode.addEdge(input.get(0), input.get(1));
		KmerPathNode.addEdge(input.get(0), input.get(2));
		KmerPathNode.addEdge(input.get(2), input.get(3));
		KmerPathNode.addEdge(input.get(3), input.get(4));
		List<KmerPathNode> result = simplify(go(k, 100, 100, input));
		assertEquals(1, result.size());
	}
	@Test
	public void should_collapse_multi_path_split_leaf_bwd() {
		int k = 4;
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(k, "ATTCC", 20, 30, false));
		input.add(KPN(k, "GGTGACCAATT", 12, 22, false));
		input.add(KPN(k, "GGATT", 19, 27, false));
		input.add(KPN(k, "TGGGA", 17, 25, false));
		KmerPathNode.addEdge(input.get(1), input.get(0));
		KmerPathNode.addEdge(input.get(2), input.get(0));
		KmerPathNode.addEdge(input.get(3), input.get(2));
		List<KmerPathNode> result = simplify(go(k, 100, 100, input));
		assertEquals(5, result.size());
	}
	@Test
	public void should_collapse_simple_path() {
		int k = 4;
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(k, "GTAC", 2, 10, false)); 
		input.add(KPN(k, "TACTAAA", 3, 11, false, 2));
		input.add(KPN(k, "TACGAAA", 3, 11, false));
		input.add(KPN(k, "AAAT", 7, 15, false));
		KmerPathNode.addEdge(input.get(0), input.get(1));
		KmerPathNode.addEdge(input.get(0), input.get(2));
		KmerPathNode.addEdge(input.get(1), input.get(3));
		KmerPathNode.addEdge(input.get(2), input.get(3));
		List<KmerPathNode> result = simplify(go(k, 100, 100, input));
		assertEquals(1, result.size());
	}
	@Test
	public void should_not_collapse_nonoverlapping() {
		int k = 4;
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(k, "GTAC", 2, 10, false)); 
		input.add(KPN(k, "TACTAAA", 3, 4, false));
		input.add(KPN(k, "TACGAAA", 5, 6, false, 2));
		KmerPathNode.addEdge(input.get(0), input.get(1));
		KmerPathNode.addEdge(input.get(0), input.get(2));
		List<KmerPathNode> result = simplify(go(k, 100, 100, input));
		assertEquals(3, result.size());
	}
	@Test
	public void should_not_collapse_nonoverlapping_with_common_path() {
		int k = 4;
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(k, "GTAC", 2, 10, false)); 
		input.add(KPN(k, "TACTAAA", 3, 4, false));
		input.add(KPN(k, "TACGAAA", 5, 6, false, 2));
		input.add(KPN(k, "TACCCTG", 3, 11, false, 2));
		KmerPathNode.addEdge(input.get(0), input.get(1));
		KmerPathNode.addEdge(input.get(0), input.get(2));
		KmerPathNode.addEdge(input.get(0), input.get(3));
		List<KmerPathNode> result = simplify(go(k, 100, 3, input));
		assertEquals(4, result.size());
	}
	@Test
	public void should_collapse_narrowed_interval() {
		int k = 4;
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(k, "GTAC", 2, 10, false)); 
		input.add(KPN(k, "TACT", 5, 6, false)); // gets removed -1 
		input.add(KPN(k, "ACTG", 4, 12, false)); // gets split into three separate nodes with the middle one removed + 1
		input.add(KPN(k, "TACCG", 3, 11, false, 2)); // gets split into three separate nodes with the middle one having added weight + 2
		KmerPathNode.addEdge(input.get(0), input.get(1));
		KmerPathNode.addEdge(input.get(1), input.get(2));
		KmerPathNode.addEdge(input.get(0), input.get(3));
		List<KmerPathNode> result = simplify(go(k, 100, 3, input));
		assertEquals(6, result.size());
	}
	@Test
	public void should_not_collapse_non_overlapping_end_paths() {
		int k = 4;
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(k, "GTAC", 1, 2, false)); 
		input.add(KPN(k, "TACT", 2, 3, false)); 
		input.add(KPN(k, "ACTG", 3, 4, false));
		input.add(KPN(k, "CTGATT", 5, 5, false));
		input.add(KPN(k, "TACC", 2, 3, false, 2)); 
		input.add(KPN(k, "ACCG", 3, 4, false, 2));
		input.add(KPN(k, "CCGTTT", 4, 4, false, 2));
		// add these so we don't have overlapping shorter paths
		input.add(KPN(k, "CTGGCGTAGA", 4, 4, false));
		input.add(KPN(k, "CCGGGCAACT", 5, 5, false));
		KmerPathNode.addEdge(input.get(0), input.get(1));
		KmerPathNode.addEdge(input.get(1), input.get(2));
		KmerPathNode.addEdge(input.get(2), input.get(3));
		KmerPathNode.addEdge(input.get(0), input.get(4));
		KmerPathNode.addEdge(input.get(4), input.get(5));
		KmerPathNode.addEdge(input.get(5), input.get(6));
		KmerPathNode.addEdge(input.get(2), input.get(7));
		KmerPathNode.addEdge(input.get(5), input.get(8));
		
		List<KmerNode> splitExpected = split(input);
		List<KmerPathNode> result = go(k, 100, 2, input);
		List<KmerNode> splitActual = split(result);
		assertEquals(splitExpected, splitActual);
	}
}
