package au.edu.wehi.idsv.debruijn.positional;

import static org.junit.Assert.assertEquals;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import org.junit.Test;

import au.edu.wehi.idsv.DiscordantReadPair;
import au.edu.wehi.idsv.SoftClipEvidence;
import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;

import com.google.common.collect.Lists;


public class BestNonReferenceContigCallerTest extends TestHelper {
	public List<ArrayDeque<KmerPathSubnode>> calls(List<KmerPathNode> input, int maxEvidenceWidth) {
		List<ArrayDeque<KmerPathSubnode>> list = new ArrayList<ArrayDeque<KmerPathSubnode>>();
		BestNonReferenceContigCaller caller = new BestNonReferenceContigCaller(input.iterator(), 100);
		for (ArrayDeque<KmerPathSubnode> contig = caller.bestContig(); contig != null; contig = caller.bestContig()) {
			list.add(contig);
		}
		return list;
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
		assertEquals(1, result.get(0).getFirst().firstKmerStartPosition());
		assertEquals(10, result.get(0).getFirst().firstKmerEndPosition());
		// alternate start position
		assertEquals(1, result.get(1).size());
		assertEquals(12, result.get(1).getFirst().firstKmerStartPosition());
		assertEquals(13, result.get(1).getFirst().firstKmerEndPosition());
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
		List<String> result = contigs(Lists.newArrayList(asKPN(4, 100, 100, 100, spanf, spanb, scf, scb, e)), 100, 4);		
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
		DiscordantReadPair rp = (DiscordantReadPair)NRRP(SES(10, 200), withSequence("GACCTCCGGAA", DP(0, 25, "10M", true, 1, 1, "10M", false)));
		ArrayList<KmerPathNode> in = Lists.newArrayList(asKPN(4, 1000, 1000, 1000, sc, rp));
		List<String> result = contigs(in, 1000, 4);
		assertEquals(3, result.size()); // SC+RP, RP starting before SC, RP starting after SC
		assertEquals("TGGTCGACCTCCGGAA", result.get(0));
		assertEquals("GACCTCCGGAA", result.get(1));
		assertEquals("GACCTCCGGAA", result.get(2));
	}
}
