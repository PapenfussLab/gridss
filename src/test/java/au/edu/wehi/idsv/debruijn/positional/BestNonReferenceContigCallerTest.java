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

import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;


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
}
