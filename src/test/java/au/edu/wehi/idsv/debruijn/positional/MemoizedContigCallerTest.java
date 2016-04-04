package au.edu.wehi.idsv.debruijn.positional;

import java.util.ArrayList;

import com.google.api.client.util.Lists;


public class MemoizedContigCallerTest extends ContigCallerTest {
	@Override
	public ContigCaller getCaller(Iterable<KmerPathNode> input, int maxEvidenceWidth) {
		ArrayList<KmerPathNode> x = Lists.newArrayList(input);
		x.sort(KmerNodeUtil.ByFirstStart);
		return new MemoizedContigCaller(x.iterator(), maxEvidenceWidth);
	}

}
