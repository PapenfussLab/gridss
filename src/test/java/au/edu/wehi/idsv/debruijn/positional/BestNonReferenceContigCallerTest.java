package au.edu.wehi.idsv.debruijn.positional;

import static org.junit.Assert.assertEquals;

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

import au.edu.wehi.idsv.DiscordantReadPair;
import au.edu.wehi.idsv.SoftClipEvidence;
import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.util.IntervalUtil;

import com.google.common.collect.Lists;


public class BestNonReferenceContigCallerTest extends TestHelper {
	public List<ArrayDeque<KmerPathSubnode>> calls(List<KmerPathNode> input, int maxEvidenceWidth) {
		BestNonReferenceContigCaller caller = new BestNonReferenceContigCaller(input.iterator(), maxEvidenceWidth);
		List<ArrayDeque<KmerPathSubnode>> contigs = caller.contigsFound();
		return contigs;
	}
	public List<String> contigs(List<KmerPathNode> input, int maxEvidenceWidth, int k) {
		Function<KmerPathSubnode, Stream<Long>> toKmers = sn -> IntStream.range(0, sn.length()).mapToObj(i -> sn.kmer(i));
		return calls(input, maxEvidenceWidth).stream().map(
					contig -> S(KmerEncodingHelper.baseCalls(contig.stream().flatMap(toKmers).collect(Collectors.toList()), k))
				).collect(Collectors.toList());
	}
	@Test
	public void should_return_best_path_first() {
		int k = 4;
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(k, "AAAA", 1, 10, false));
		input.add(KPN(k, "AAAT", 2, 13, false));
		KmerPathNode.addEdge(input.get(0), input.get(1));
		List<ArrayDeque<KmerPathSubnode>> result = calls(input, 100);
		assertEquals(2, result.size());
		assertEquals(2, result.get(0).size());
		assertEquals(1, result.get(0).getFirst().firstStart());
		assertEquals(10, result.get(0).getFirst().firstEnd());
		// alternate start position
		assertEquals(1, result.get(1).size());
		assertEquals(12, result.get(1).getFirst().firstStart());
		assertEquals(13, result.get(1).getFirst().firstEnd());
	}
	@Test
	public void should_start_path_with_no_predecessor() {
		int k = 4;
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(k, "AAAA", 1, 10, false));
		input.add(KPN(k, "AAAT", 2, 11, false));
		KmerPathNode.addEdge(input.get(0), input.get(1));
		List<String> result = contigs(input, 100, k);
		assertEquals(1, result.size());
		assertEquals("AAAAT", result.get(0));
	}
	@Test
	public void should_not_start_in_reference_kmer() {
		int k = 4;
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(k, "AAAA", 1, 10, true));
		List<String> result = contigs(input, 100, k);
		assertEquals(0, result.size());
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
		List<String> result = contigs(Lists.newArrayList(asKPN(4, 100, spanf, spanb, scf, scb, e)), 100, 4);		
		assertEquals(4, result.size());
		// only returns the unanchored kmers
		assertEquals("GTTAATTG", result.get(0));
		assertEquals("TGGTCGACC", result.get(1));
		assertEquals("GTCAGT", result.get(2));
		assertEquals("GACCTCTACT", result.get(3));
	}
	@Test
	public void should_assemble_overlapping_sc_rp() {
		SoftClipEvidence sc = SCE(FWD, withSequence("ACGTGGTCGACC", Read(0, 50, "6M6S")));
		DiscordantReadPair rp = (DiscordantReadPair)NRRP(SES(10, 200), withSequence("GACCTCCGGAA", DP(0, 25, "11M", true, 1, 1, "11M", false)));
		ArrayList<KmerPathNode> in = Lists.newArrayList(asKPN(4, 1000, sc, rp));
		List<String> result = contigs(in, 1000, 4);
		assertEquals(3, result.size()); // SC+RP, RP starting before SC, RP starting after SC
		assertEquals("TGGTCGACCTCCGGAA", result.get(0));
		assertEquals("GACCTCCGGAA", result.get(1));
		assertEquals("GACCTCCGGAA", result.get(2));
	}
	@Test
	public void should_memoize_paths() {
		List<KmerPathNode> list = new ArrayList<KmerPathNode>();
		for (int len = 1; len <= 4; len++) {
			for (int start = 1; start <= 4; start++) {
				for (int end = start; end <= 5; end++) {
					list.add(KPN(4,StringUtils.repeat("A", 3 + len), start, end, false, start));
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
		list.sort(KmerNodeUtil.ByFirstStart);
		BestNonReferenceContigCaller caller = new BestNonReferenceContigCaller(list.iterator(), 3);
		// 1 -> 2 -> 3 -> 4 -> 6[4-6]
		ArrayDeque<KmerPathSubnode> contig = caller.bestContig();
		assertEquals(list.size(), caller.tracking_memoizedNodeCount());
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
			list.add(KPN(4, StringUtils.repeat("A", 4 + rng.nextInt(10)), start, end, false, weight));
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
		BestNonReferenceContigCaller caller = new BestNonReferenceContigCaller(list.iterator(), 3);
		caller.bestContig();
		assertEquals(list.size(), caller.tracking_memoizedNodeCount());
		assertEquals(0, caller.tracking_frontierSize());
	}
	@Test
	public void should_traverse_subnodes_through_node_at_multiple_positions() {
		int k = 4;
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(k, "ATAT", 1, 11, false));
		input.add(KPN(k, "TATA", 2, 10, false));
		// best contig:
		// 1 3 5 7 9  11
		//  2 4 6 8 10
		// second best:
		// 2 4 6 8 10
		//  3 5 7 9  
		BestNonReferenceContigCaller caller = new BestNonReferenceContigCaller(input.iterator(), 3);
		ArrayDeque<KmerPathSubnode> best = caller.bestContig();
		assertEquals(11, best.size());
	}
	@Test
	public void should_self_traverse_node() {
		int k = 4;
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(k, "AAAA", 1, 11, false));
		// best contig:
		// 1 3 5 7 9  11
		//  2 4 6 8 10
		BestNonReferenceContigCaller caller = new BestNonReferenceContigCaller(input.iterator(), 3);
		ArrayDeque<KmerPathSubnode> best = caller.bestContig();
		assertEquals(11, best.size());
	}
	@Test
	public void subnode_loops_should_not_cause_stack_overflow() {
		int k = 3;
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(k, "ACT", 1, 100000, false));
		input.add(KPN(k, "CTA", 1, 100000, false));
		input.add(KPN(k, "TAC", 1, 100000, false));
		BestNonReferenceContigCaller caller = new BestNonReferenceContigCaller(input.iterator(), 1000000);
		ArrayDeque<KmerPathSubnode> best = caller.bestContig();
		assertEquals(100000, best.size());
	}
}




