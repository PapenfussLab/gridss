package au.edu.wehi.idsv.debruijn.positional;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import org.apache.commons.lang3.StringUtils;
import org.junit.Test;

import com.google.common.collect.Lists;

import au.edu.wehi.idsv.DiscordantReadPair;
import au.edu.wehi.idsv.SoftClipEvidence;
import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.util.IntervalUtil;

public abstract class ContigCallerTest extends TestHelper {
	public abstract ContigCaller getCaller(Iterable<KmerPathNode> input, int maxEvidenceWidth);
	public static String S(ArrayDeque<KmerPathSubnode> contig, int k) {
		Function<KmerPathSubnode, Stream<Long>> toKmers = sn -> IntStream.range(0, sn.length()).mapToObj(i -> sn.kmer(i));
		return S(KmerEncodingHelper.baseCalls(contig.stream().flatMap(toKmers).collect(Collectors.toList()), k));
	}
	@Test
	public void should_memoize_paths() {
		List<KmerPathNode> list = new ArrayList<KmerPathNode>();
		for (int len = 1; len <= 4; len++) {
			for (int start = 1; start <= 4; start++) {
				for (int end = start; end <= 5; end++) {
					// ensure kmers are unique for each starting position
					long kmer = 6 * ((5 * len) + start) + end;
					String startingKmer = KmerEncodingHelper.toString(10, kmer);
					list.add(KPN(10, startingKmer + StringUtils.repeat("A", len - 1), start, end, false, start));
				}
			}
		}
		//list.add(KPN(4, "AAAA", 1, 100, false, 1));
		for (KmerPathNode prev : list) {
			for (KmerPathNode next : list) {
				if (IntervalUtil.overlapsClosed(prev.lastStart() + 1, prev.lastEnd() + 1, next.firstStart(), next.firstEnd())) {
					KmerPathNode.addEdge(prev, next);
				}
			}
		}
		ContigCaller caller = getCaller(list, 10);
		// highest scoring path goes through the most nodes:
		// [1-1][1-*] -> [2-2][2-*] -> [3-3][3-*] -> [4-4][4-*] -> [5-5][4-*] 
		ArrayDeque<KmerPathSubnode> contig = caller.bestContig(Integer.MAX_VALUE);
		assertEquals(list.size(), caller.tracking_memoizedNodeCount());
		assertEquals(0, caller.tracking_frontierSize());
		assertEquals(5, contig.size());
	}
	@Test
	public void should_memoize_paths_random_path_test() {
		// assertion testing
		Random rng = new Random(0);
		List<KmerPathNode> list = new ArrayList<KmerPathNode>();
		for (int i = 0; i < 32; i++) {
			int start = 1 + rng.nextInt(100);
			int end = start + rng.nextInt(100);
			int weight = 1 + rng.nextInt(100);
			// ensure kmers are unique for each starting position
			long kmer = 100 * ((100 * weight) + start) + end;
			String startingKmer = KmerEncodingHelper.toString(32, kmer);
			
			list.add(KPN(32, startingKmer + StringUtils.repeat("A", rng.nextInt(10)), start, end, false, weight));
		}
		//list.add(KPN(4, "AAAA", 1, 100, false, 1));
		for (KmerPathNode prev : list) {
			for (KmerPathNode next : list) {
				if (IntervalUtil.overlapsClosed(prev.lastStart() + 1, prev.lastEnd() + 1, next.firstStart(), next.firstEnd())) {
					KmerPathNode.addEdge(prev, next);
				}
			}
		}
		list.sort(KmerNodeUtil.ByFirstStart);
		ContigCaller caller = getCaller(list, 3);
		caller.bestContig(Integer.MAX_VALUE);
		assertEquals(list.size(), caller.tracking_memoizedNodeCount());
		assertEquals(0, caller.tracking_frontierSize());
	}
	@Test
	public void should_traverse_subnodes_through_node_at_multiple_positions() {
		int k = 4;
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(k, "ATAT", 1, 11, false));
		input.add(KPN(k, "TATA", 2, 10, false));
		KmerPathNode.addEdge(input.get(0), input.get(1));
		KmerPathNode.addEdge(input.get(1), input.get(0));
		// best contig:
		// 1 3 5 7 9  11
		//  2 4 6 8 10
		// second best:
		// 2 4 6 8 10
		//  3 5 7 9  
		ContigCaller caller = getCaller(input, 3);
		ArrayDeque<KmerPathSubnode> best = caller.bestContig(Integer.MAX_VALUE);
		assertEquals(11, best.size());
	}
	@Test
	public void should_self_traverse_node() {
		int k = 4;
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(k, "AAAA", 1, 11, false));
		KmerPathNode.addEdge(input.get(0), input.get(0));
		// best contig:
		// 1 3 5 7 9  11
		//  2 4 6 8 10
		ContigCaller caller = getCaller(input, 3);
		ArrayDeque<KmerPathSubnode> best = caller.bestContig(Integer.MAX_VALUE);
		assertEquals(11, best.size());
	}
	@Test
	public void subnode_loops_should_not_cause_stack_overflow() {
		int length = 1000;
		int k = 3;
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(k, "ACT", 1, length, false));
		input.add(KPN(k, "CTA", 1, length, false));
		input.add(KPN(k, "TAC", 1, length, false));
		KmerPathNode.addEdge(input.get(0), input.get(1));
		KmerPathNode.addEdge(input.get(1), input.get(2));
		KmerPathNode.addEdge(input.get(2), input.get(0));
		ContigCaller caller = getCaller(input, length);
		ArrayDeque<KmerPathSubnode> best = caller.bestContig(Integer.MAX_VALUE);
		assertEquals(length, best.size());
	}
	@Test
	public void should_start_path_with_no_predecessor() {
		int k = 4;
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(k, "AAAA", 1, 10, false));
		input.add(KPN(k, "AAAT", 2, 11, false));
		KmerPathNode.addEdge(input.get(0), input.get(1));
		ContigCaller caller = getCaller(input, 100000);
		String result = S(caller.bestContig(Integer.MAX_VALUE), k);
		assertEquals("AAAAT", result);
	}
	@Test
	public void should_not_start_in_reference_kmer() {
		int k = 4;
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(k, "AAAA", 1, 10, true));
		ArrayDeque<KmerPathSubnode> result = getCaller(input, 100000).bestContig(Integer.MAX_VALUE);
		assertNull(result);
	}
	@Test
	public void should_return_best_path() {
		int k = 4;
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(k, "AAAA", 1, 10, false));
		input.add(KPN(k, "AAAT", 2, 13, false));
		KmerPathNode.addEdge(input.get(0), input.get(1));
		ArrayDeque<KmerPathSubnode> result = getCaller(input, 10).bestContig(Integer.MAX_VALUE);
		assertEquals(2, result.size());
		assertEquals(1, result.getFirst().firstStart());
		assertEquals(10, result.getFirst().firstEnd());
		// alternate start position
		//assertEquals(1, result.size());
		//assertEquals(12, result.getFirst().firstStart());
		//assertEquals(13, result.getFirst().firstEnd());
	}
	@Test
	public void should_assemble_overlapping_sc_rp() {
		SoftClipEvidence sc = SCE(FWD, withSequence("ACGTGGTCGACC", Read(0, 50, "6M6S")));
		DiscordantReadPair rp = (DiscordantReadPair)NRRP(SES(10, 200), withSequence("GACCTCCGGAA", DP(0, 25, "11M", true, 1, 1, "11M", false)));
		ArrayList<KmerPathNode> in = Lists.newArrayList(asKPN(4, 1000, sc, rp));
		String result = S(getCaller(in, 1000).bestContig(Integer.MAX_VALUE), 4);
		//assertEquals(3, result.size()); // SC+RP, RP starting before SC, RP starting after SC
		assertEquals("TGGTCGACCTCCGGAA", result);
		//assertEquals("GACCTCCGGAA", result.get(1));
		//assertEquals("GACCTCCGGAA", result.get(2));
	}
	@Test
	public void should_preference_anchored_paths() {
		// both sides anchored
		SoftClipEvidence spanf = SCE(FWD, withSequence("TGTTAATTGT", Read(0, 1, "4M6S")));
		SoftClipEvidence spanb = SCE(BWD, withSequence("TGTTAATTGT", Read(0, 7, "6S4M")));
		// start anchored
		SoftClipEvidence scf = SCE(FWD, withSequence("ACGTGGTCGACC", Read(0, 5, "6M6S")));
		// end anchored (but less support)
		SoftClipEvidence scb = SCE(BWD, withSequence("GTCAGTC", Read(0, 5, "3S4M")));
		// unanchored
		DiscordantReadPair e = (DiscordantReadPair)NRRP(SES(20, 30), withSequence("GACCTCTACT", DP(0, 25, "10M", true, 1, 1, "10M", false)));
		
		String result = S(getCaller(Lists.newArrayList(asKPN(4, 100, spanf, spanb, scf, scb, e)), 100).bestContig(Integer.MAX_VALUE), 4);
		//assertEquals(4, result.size());
		// only returns the unanchored kmers
		assertEquals("GTTAATTG", result);
		//assertEquals("TGGTCGACC", result.get(1));
		//assertEquals("GTCAGT", result.get(2));
		//assertEquals("GACCTCTACT", result.get(3));
	}
	@Test
	public void should_not_call_suboptimal_contig_with_successor_within_evidence_distance() {
		List<KmerPathNode> list = new ArrayList<KmerPathNode>();
		list.add(KPN(4, "TTTT", 1, 1, false, 70));
		for (int i = 10; i < 110; i++) {
			list.add(KPN(4, "AAAA", i, i, false, 1));
		}
		for (int i = 2; i < list.size(); i++) {
			KmerPathNode.addEdge(list.get(i - 1), list.get(i));
		}
		//          1         2         3
		// 123456789012345678901234567890
		// X        XXXXXXXXXXXXXXXXXXXX
		//
		ArrayDeque<KmerPathSubnode> result = getCaller(list, 10).bestContig(Integer.MAX_VALUE);
		assertEquals(100, result.size());
	}
	@Test
	public void should_favor_early_calls() {
		List<KmerPathNode> list = new ArrayList<KmerPathNode>();
		list.add(KPN(4, "TTTT", 2, 2, false, 100));
		list.add(KPN(4, "TTTT", 3, 3, false, 100));
		ContigCaller caller = getCaller(list, 2);
		caller.bestContig(Integer.MAX_VALUE);
		caller.add(KPN(4, "TTTT", 4, 4, false, 100));
		caller.add(KPN(4, "TTTT", 1, 1, false, 100));
		assertEquals(1, caller.bestContig(Integer.MAX_VALUE).getFirst().firstStart());
	}
	@Test
	public void should_call_first_optimal_contig() {
		List<KmerPathNode> list = new ArrayList<KmerPathNode>();
		list.add(KPN(4, "TTTT", 1, 1, false, 20));
		for (int i = 10; i < 30; i++) {
			list.add(KPN(4, "AAAA", i, i, false, 1));
		}
		for (int i = 2; i < list.size(); i++) {
			KmerPathNode.addEdge(list.get(i - 1), list.get(i));
		}
		//          1         2         3
		// 123456789012345678901234567890
		// X        XXXXXXXXXXXXXXXXXXXX
		//
		ArrayDeque<KmerPathSubnode> result = getCaller(list, 2).bestContig(Integer.MAX_VALUE);
		assertEquals(1, result.size());
	}
	@Test
	public void should_not_call_contig_if_unclear_if_optimal() {
		List<KmerPathNode> list = new ArrayList<KmerPathNode>();
		list.add(KPN(4, "TTTT", 1, 1, false, 10));
		for (int i = 10; i < 30; i++) {
			list.add(KPN(4, "AAAA", i, i, false, 1));
		}
		for (int i = 2; i < list.size(); i++) {
			KmerPathNode.addEdge(list.get(i - 1), list.get(i));
		}
		//          1         2         3
		// 123456789012345678901234567890
		// X        XXXXXXXXXXXXXXXXXXXX
		//
		ArrayDeque<KmerPathSubnode> result = getCaller(list, 2).bestContig(28);
		assertNull(result);
	}
	@Test
	public void should_not_call_reference_paths() {
		int k = 4;
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(k, "AAAA", 1, 1, true));
		input.add(KPN(k, "AAAA", 2, 2, true));
		ArrayDeque<KmerPathSubnode> result = getCaller(input, 10).bestContig(Integer.MAX_VALUE);
		assertNull(result);
	}
	@Test
	public void reference_nodes_should_be_terminal() {
		int k = 4;
		List<KmerPathNode> list = new ArrayList<KmerPathNode>();
		for (int i = 1; i < 16; i++) {
			list.add(KPN(k, "AAAA", i, i, i % 2 == 0));
		}
		for (int i = 1; i < list.size(); i++) {
			KmerPathNode.addEdge(list.get(i - 1), list.get(i));
		}
		ArrayDeque<KmerPathSubnode> result = getCaller(list, 10).bestContig(Integer.MAX_VALUE);
		assertEquals(1, result.size());
		assertTrue(result.getFirst().prev().get(0).isReference());
		assertTrue(result.getLast().next().get(0).isReference());
	}
	@Test
	public void should_not_include_anchor_weight() {
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(4, "TTTT", 1, 1, false, 2));
		input.add(KPN(4, "TTTT", 2, 2, true, 1));
		input.add(KPN(4, "AAAA", 3, 3, false, 1));
		input.add(KPN(4, "AAAA", 4, 4, true, 10));
		KmerPathNode.addEdge(input.get(0), input.get(1));
		KmerPathNode.addEdge(input.get(2), input.get(3));
		ArrayDeque<KmerPathSubnode> result = getCaller(input, 10).bestContig(Integer.MAX_VALUE);
		// should take the path that has higher non-reference weight
		assertEquals(1, result.getFirst().firstStart());
	}
	@Test
	public void should_call_path_joining_reference_path() {
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(4, "AAAA", 1, 1, false, 1));
		input.add(KPN(4, "TAAA", 1, 1, true, 1));
		input.add(KPN(4, "AAAA", 2, 2, true, 1));
		KmerPathNode.addEdge(input.get(0), input.get(2));
		KmerPathNode.addEdge(input.get(1), input.get(2));
		ArrayDeque<KmerPathSubnode> result = getCaller(input, 10).bestContig(Integer.MAX_VALUE);
		// should take the path that has higher non-reference weight
		assertEquals(1, result.getFirst().firstStart());
	}
}
