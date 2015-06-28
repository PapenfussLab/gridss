package au.edu.wehi.idsv.debruijn.positional;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import org.apache.commons.lang3.ArrayUtils;
import org.junit.Test;

import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.DiscordantReadPair;
import au.edu.wehi.idsv.SoftClipEvidence;
import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.util.IntervalUtil;

import com.google.common.collect.Lists;


public class PathNodeIteratorTest extends TestHelper {
	/**
	 * Checks that all edges are included in the given graph
	 * @param list
	 * @param k
	 */
	public static void assertCompleteGraph(List<KmerPathNode> list, final int k) {
		for (KmerPathNode pn : list) {
			long[] nextKmers = KmerEncodingHelper.nextStates(k, pn.lastKmer());
			List<KmerPathNode> nextList = list.stream()
					.filter(n -> ArrayUtils.contains(nextKmers, n.kmer(0)) && IntervalUtil.overlapsClosed(n.startPosition(0), n.endPosition(0), pn.lastStart() + 1, pn.lastEnd() + 1))
					.collect(Collectors.toList());
			assertEquals(nextList.size(), pn.next().size());
			assertTrue(KmerNodeUtil.ByFirstStart.isOrdered(pn.next()));
			nextList.stream().forEach(n -> assertTrue(pn.next().contains(n)));
			
			long[] prevKmers = KmerEncodingHelper.prevStates(k, pn.kmer(0));
			List<KmerPathNode> prevList = list.stream()
					.filter(n -> ArrayUtils.contains(prevKmers, n.lastKmer()) && IntervalUtil.overlapsClosed(n.lastStart(), n.lastEnd(), pn.startPosition(0) - 1, pn.endPosition(0) - 1))
					.collect(Collectors.toList());
			assertEquals(prevList.size(), pn.prev().size());
			assertTrue(KmerNodeUtil.ByLastStart.isOrdered(pn.prev()));
			prevList.stream().forEach(n -> assertTrue(pn.prev().contains(n)));
		}
	}
	private void assertContains(List<KmerPathNode> list, int k, String baseSequence, int firstStart, int firstEnd, boolean reference, int weight) {
		for (KmerPathNode node : list) {
			if (matches(node, k, baseSequence, firstStart, firstEnd, reference, weight)) {
				return;
			}
		}
		fail(String.format("Missing (%s, %d, %d, %s, %d)", baseSequence, firstStart, firstEnd, reference, weight));
	}
	private boolean matches(KmerPathNode node, int k, String baseSequence, int firstStart, int firstEnd, boolean reference, int weight) {
		if (firstStart != node.startPosition(0)) return false;
		if (firstEnd != node.endPosition(0)) return false;
		if (weight != node.weight()) return false;
		if (node.isReference() != reference) return false;
		List<Long> list = toKmers(baseSequence, k);
		if (list.size() != node.length()) return false;
		for (int i = 0; i < list.size(); i++) {
			if ((long)list.get(i) != node.kmer(i)) return false;
		}
		return true;
	}
	private List<Long> toKmers(String baseSequence, int k) {
		List<Long> list = new ArrayList<Long>();
		for (int i = 0; i + k <= baseSequence.length(); i++) {
			list.add(K(baseSequence.substring(i, i + k)));
		}
		return list;
	}
	/**
	 * Create a data set that tiles
	 */
	public static List<KmerNode> tiled(int width, int tilesWide, int length, int k) {
		List<KmerNode> list = new ArrayList<KmerNode>();
		for (int lengthOffset = 0; lengthOffset < length; lengthOffset++) {
			for (int i = 0; i < tilesWide; i++) {
				int start = i * width;
				int end = start + width - 1;			
				list.add(new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, Arrays.copyOfRange(RANDOM, 0+lengthOffset, k+lengthOffset)), start+lengthOffset, end+lengthOffset, true, 1));
			}
		}
		list.sort(KmerNodeUtil.ByLastStart);
		// sanity checks
		assertTrue(list.stream().allMatch(n -> n.width() == width));
		assertTrue(list.size() == list.stream().distinct().count());
		assertDisjoint(list);
		return list;
	}
	@Test
	public void should_merge_successive_nodes() {
		List<ImmutableKmerNode> in = new ArrayList<ImmutableKmerNode>();
		in.add(new ImmutableKmerNode(K("GAAA"), 1, 2, true, 1));
		in.add(new ImmutableKmerNode(K("AAAT"), 2, 3, true, 1));
		List<KmerPathNode> out = Lists.newArrayList(new PathNodeIterator(in.iterator(), 10, 4));
		assertCompleteGraph(out, 4);
		assertEquals(1, out.size());
		assertContains(out, 4, "GAAAT", 1, 2, true, 2);
		assertSameNodes(in, out);
	}
	public void assertCount(int expected, ImmutableKmerNode... in) {
		List<ImmutableKmerNode> list = Lists.newArrayList(in);
		Collections.sort(list, KmerNodeUtil.ByLastStart);
		List<KmerPathNode> out = Lists.newArrayList(new PathNodeIterator(list.iterator(), 10, 4));
		assertCompleteGraph(out, 4);
		assertEquals(expected, out.size());
		assertSameNodes(list, out);
	}
	@Test
	public void should_only_merge_fully_matching_nodes() {
		assertCount(1,
			new ImmutableKmerNode(K("GAAA"), 1, 2, true, 1),
			new ImmutableKmerNode(K("AAAT"), 2, 3, true, 1));
		assertCount(2,
				new ImmutableKmerNode(K("GAAA"), 1, 2, true, 1),
				new ImmutableKmerNode(K("AAAT"), 2, 3, false, 1));
		assertCount(2,
				new ImmutableKmerNode(K("GAAA"), 1, 3, true, 1),
				new ImmutableKmerNode(K("AAAT"), 2, 3, true, 1));
		assertCount(2,
				new ImmutableKmerNode(K("GAAA"), 1, 2, true, 1),
				new ImmutableKmerNode(K("AAAT"), 2, 4, true, 1));
		assertCount(1,
				new ImmutableKmerNode(K("GAAA"), 1, 2, true, 1),
				new ImmutableKmerNode(K("AAAT"), 2, 3, true, 2));
	}
	@Test
	public void should_handle_self_intersection() {
		assertCount(1, new ImmutableKmerNode(K("AAAA"), 1, 5, true, 1));
	}
	@Test
	public void should_split_path_upon_fork() {
		List<ImmutableKmerNode> in = new ArrayList<ImmutableKmerNode>();
		in.add(new ImmutableKmerNode(K("TAAA"), 0, 0, true, 1));
		in.add(new ImmutableKmerNode(K("AAAA"), 1, 1, true, 1));
		in.add(new ImmutableKmerNode(K("AAAA"), 2, 2, true, 1));
		in.add(new ImmutableKmerNode(K("AAAT"), 2, 2, true, 1));
		in.add(new ImmutableKmerNode(K("AAAG"), 2, 2, true, 1));
		in.add(new ImmutableKmerNode(K("AAAC"), 2, 2, true, 1));
		List<KmerPathNode> out = Lists.newArrayList(new PathNodeIterator(in.iterator(), 10, 4));
		assertEquals(5, out.size());
		assertContains(out, 4, "TAAAA", 0, 0, true, 2);
		assertContains(out, 4, "AAAA", 2, 2, true, 1);
		assertContains(out, 4, "AAAC", 2, 2, true, 1);
		assertContains(out, 4, "AAAG", 2, 2, true, 1);
		assertContains(out, 4, "AAAT", 2, 2, true, 1);
		assertCompleteGraph(out, 4);
	}
	@Test
	public void should_concatenate_to_path() {
		List<ImmutableKmerNode> in = new ArrayList<ImmutableKmerNode>();
		in.add(new ImmutableKmerNode(K("TAAA"), 1, 3, true, 1));
		in.add(new ImmutableKmerNode(K("AAAC"), 2, 4, true, 2));
		in.add(new ImmutableKmerNode(K("AACT"), 3, 5, true, 3));
		in.add(new ImmutableKmerNode(K("ACTT"), 4, 6, true, 4));
		in.add(new ImmutableKmerNode(K("CTTG"), 5, 7, true, 5));
		in.add(new ImmutableKmerNode(K("TTGG"), 6, 8, true, 6));
		List<KmerPathNode> out = Lists.newArrayList(new PathNodeIterator(in.iterator(), 10, 4));
		assertCompleteGraph(out, 4);
		assertEquals(1, out.size());
		assertContains(out, 4, "TAAACTTGG", 1, 3, true, 1+2+3+4+5+6);
	}
	@Test
	public void should_limit_path_length() {
		List<ImmutableKmerNode> in = new ArrayList<ImmutableKmerNode>();
		in.add(new ImmutableKmerNode(K("TAAA"), 1, 3, true, 1));
		in.add(new ImmutableKmerNode(K("AAAC"), 2, 4, true, 2));
		in.add(new ImmutableKmerNode(K("AACT"), 3, 5, true, 3));
		in.add(new ImmutableKmerNode(K("ACTT"), 4, 6, true, 4));
		in.add(new ImmutableKmerNode(K("CTTG"), 5, 7, true, 5));
		in.add(new ImmutableKmerNode(K("TTGG"), 6, 8, true, 6));
		List<KmerPathNode> out = Lists.newArrayList(new PathNodeIterator(in.iterator(), 3, 4));
		assertCompleteGraph(out, 4);
		assertEquals(2, out.size());
		assertContains(out, 4, "TAAACT", 1, 3, true, 1+2+3);
		assertContains(out, 4, "ACTTGG", 4, 6, true, 4+5+6);
	}
	@Test
	public void should_not_collapse_reference_with_non_reference() {
		List<ImmutableKmerNode> in = new ArrayList<ImmutableKmerNode>();
		in.add(new ImmutableKmerNode(K("TAAA"), 1, 3, true, 1));
		in.add(new ImmutableKmerNode(K("AAAC"), 2, 4, false, 2));
		List<KmerPathNode> out = Lists.newArrayList(new PathNodeIterator(in.iterator(), 3, 4));
		assertCompleteGraph(out, 4);
		assertEquals(2, out.size());
	}
	@Test
	public void should_output_constant_path_node_dimensions_for_tiled_input() {
		int k = 25;
		for (int tileWidth = 1; tileWidth < 4; tileWidth++) {
			for (int nTilesWide = 1; nTilesWide < 4; nTilesWide++) {
				for (int tileLength = 1; tileLength < 4; tileLength++) {
					for (int nTilesLong = 1; nTilesLong < 4; nTilesLong++) {
						List<KmerNode> tileInput = tiled(tileWidth, nTilesWide, tileLength * nTilesLong, k);
						ArrayList<KmerPathNode> result = Lists.newArrayList(new PathNodeIterator(tileInput.iterator(), tileLength, k));
						assertEquals(nTilesWide * nTilesLong, result.size());
						for (KmerPathNode pn : result) {
							assertTrue(pn.isValid());
							pn.sanityCheck(k, tileWidth, tileLength);
							assertTrue(pn.weight() > 0);
							assertEquals(tileLength, pn.length());
							assertEquals(tileWidth, pn.width());
						}
						assertCompleteGraph(result, k);
					}
				}
			}
		}
	}
	@Test
	public void link_adjacent() { 
		List<ImmutableKmerNode> in = new ArrayList<ImmutableKmerNode>();
		in.add(new ImmutableKmerNode(K("TAAA"), 5, 5, false, 1));
		in.add(new ImmutableKmerNode(K("AAAC"), 6, 6, false, 2));
		in.add(new ImmutableKmerNode(K("AACT"), 1, 10, false, 2));
		in.sort(KmerNodeUtil.ByLastStart);
		List<KmerPathNode> paths = Lists.newArrayList(new PathNodeIterator(in.iterator(), 100, 4));
		assertCompleteGraph(paths, 4);
	}
	@Test
	public void should_assemble_adjacent_sc_rp() {
		SoftClipEvidence sc = SCE(FWD, withSequence("GTGGTCGACC", Read(0, 50, "4M6S")));
		DiscordantReadPair rp = (DiscordantReadPair)NRRP(SES(10, 200), withSequence("GACCT", DP(0, 25, "5M", true, 1, 1, "5M", false)));
		List<KmerPathNode> result = asCheckedKPN(5, 1000, sc, rp);
		assertCompleteGraph(result, 5);
		assertEquals(2, result.size());
		assertEquals(1, result.stream().mapToInt(n -> n.next().size()).sum());
		assertEquals(1, result.stream().mapToInt(n -> n.prev().size()).sum());
	}
	@Test
	public void should_overlapping_adjacent_sc_rp() {
		SoftClipEvidence sc = SCE(FWD, withSequence("ACGTGGTCGACC", Read(0, 50, "6M6S")));
		DiscordantReadPair rp = (DiscordantReadPair)NRRP(SES(13, 200), withSequence("TCGACCTCCGGAA", DP(0, 25, "13M", true, 1, 1, "13M", false)));
		List<KmerPathNode> result = asCheckedKPN(5, 1000, sc, rp);
		assertCompleteGraph(result, 5);
		assertEquals(5, result.size());
		assertEquals(4, result.stream().mapToInt(n -> n.next().size()).sum());
		assertEquals(4, result.stream().mapToInt(n -> n.prev().size()).sum());
	}
	@Test
	public void should_handle_repeated_kmers() {
		SoftClipEvidence sce = SCE(FWD, withSequence("TTTTTTTTTTTTTTTT", Read(0, 5, "6M10S")));
		DirectedEvidence rp = NRRP(SES(10, 200), withSequence("CCCCCCCCCC", DP(0, 25, "10M", true, 1, 1, "10M", false)));
		
		List<KmerPathNode> result = asCheckedKPN(4, 100, sce, rp);
		assertCompleteGraph(result, 4);
		// SC: ref repeat; nonref repeat
		// RP: single width starting kmer <-> repeat cycle <-> single width ending kmer 
		assertEquals(5, result.size());
	}
}
