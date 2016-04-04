package au.edu.wehi.idsv.debruijn.positional;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import com.google.api.client.util.Lists;

import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;


public class BestNonReferenceContigCallerTest extends ContigCallerTest {
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
	@Override
	public ContigCaller getCaller(Iterable<KmerPathNode> input, int maxEvidenceWidth) {
		ArrayList<KmerPathNode> x = Lists.newArrayList(input);
		x.sort(KmerNodeUtil.ByFirstStart);
		return new BestNonReferenceContigCaller(x.iterator(), maxEvidenceWidth);
	}
}




