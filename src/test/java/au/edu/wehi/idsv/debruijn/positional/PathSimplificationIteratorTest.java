package au.edu.wehi.idsv.debruijn.positional;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;

import com.google.common.collect.Lists;


public class PathSimplificationIteratorTest extends TestHelper {
	@Test
	public void should_not_collapse_past_limits() {
		int k = 4;
		int tiles = 2;
		for (int inputWidth = 1; inputWidth <= 4; inputWidth++) {
			for (int inputLength = 1; inputLength <= 4; inputLength++) {
				for (int maxSupportWidth = inputWidth; maxSupportWidth <= 16; maxSupportWidth++) {
					for (int maxPathLength = inputLength; maxPathLength <= 16; maxPathLength++) {
						// tile evidence such that the entire input could be compressed to a single
						// node if the limits did not exist
						List<KmerNode> in = PathNodeIteratorTest.tiled(inputWidth, tiles, maxPathLength * tiles, k);
						in.sort(KmerNodeUtil.ByStartPosition);
						int weightIn = totalWeight(in);
						List<KmerPathNode> inpn = Lists.newArrayList(new PathNodeIterator(in.iterator(), inputLength, k));
						assertEquals(weightIn, totalWeight(inpn)); // precondition
						PathSimplificationIterator psi = new PathSimplificationIterator(inpn.iterator(), maxPathLength, maxSupportWidth);
						List<KmerPathNode> list = Lists.newArrayList(psi);
						assertSameNodes(in, list);
						for (KmerPathNode pn : list) {
							pn.sanityCheck(k, maxSupportWidth, maxPathLength);
							assertTrue(pn.width() <= maxSupportWidth);
							assertTrue(pn.length() <= maxPathLength);
						}
						assertEquals(weightIn, totalWeight(list));
					}
				}
			}
		}
	}
	@Test
	public void should_match_reference_when_collapsing() {
		List<KmerPathNode> in = new ArrayList<KmerPathNode>();
		in.add(new KmerPathNode(K("TAAA"), 1, 1, true, 1));
		in.add(new KmerPathNode(K("AAAT"), 2, 2, false, 1));
		KmerPathNode.addEdge(in.get(0), in.get(1));
		in.sort(KmerNodeUtil.ByFirstKmerStartPosition);
		List<KmerPathNode> list = Lists.newArrayList(new PathSimplificationIterator(in.iterator(), 64, 64));
		assertEquals(2, list.size());
		assertEquals(totalWeight(in), totalWeight(list));
	}
	@Test
	public void should_collapse_adjacent() {
		List<KmerPathNode> in = new ArrayList<KmerPathNode>();
		in.add(new KmerPathNode(K("TAAA"), 1, 1, true, 1));
		in.add(new KmerPathNode(K("TAAA"), 2, 3, true, 1));
		in.add(new KmerPathNode(K("TAAA"), 4, 4, true, 1));
		in.add(new KmerPathNode(K("TAAA"), 5, 10, true, 1));
		in.sort(KmerNodeUtil.ByFirstKmerStartPosition);
		
		int weightIn = totalWeight(in);
		List<KmerPathNode> list = Lists.newArrayList(new PathSimplificationIterator(in.iterator(), 64, 64));
		assertEquals(1, list.size());
		assertEquals(weightIn, totalWeight(list));
	}
	@Test
	public void should_collapse_consecutive() {
		List<KmerPathNode> in = new ArrayList<KmerPathNode>();
		in.add(new KmerPathNode(K("TAAA"), 1, 10, true, 1));
		in.add(new KmerPathNode(K("AAAT"), 2, 11, true, 1));
		in.get(in.size() - 1).append(new ImmutableKmerNode(K("AATC"), 3, 12, true, 5));
		in.add(new KmerPathNode(K("ATCC"), 4, 13, true, 6));
		in.sort(KmerNodeUtil.ByFirstKmerStartPosition);
		
		KmerPathNode.addEdge(in.get(0), in.get(1));
		KmerPathNode.addEdge(in.get(1), in.get(2));
		int weightIn = totalWeight(in);
		List<KmerPathNode> list = Lists.newArrayList(new PathSimplificationIterator(in.iterator(), 64, 64));
		assertEquals(1, list.size());
		assertEquals(10, list.get(0).width());
		assertEquals(4, list.get(0).length());
		assertEquals(1+1+5+6, list.get(0).weight());
		assertEquals(weightIn, totalWeight(list));
	}
	@Test
	public void should_chain_collapse() {
		List<KmerPathNode> in = new ArrayList<KmerPathNode>();
		in.add(new KmerPathNode(K("TAAA"), 1, 2, true, 4));
		in.add(new KmerPathNode(K("AAAT"), 2, 3, true, 5));
		in.add(new KmerPathNode(K("TAAA"), 3, 4, true, 4));
		in.add(new KmerPathNode(K("AAAT"), 4, 5, true, 5));
		KmerPathNode.addEdge(in.get(0), in.get(1));
		KmerPathNode.addEdge(in.get(2), in.get(3));
		in.sort(KmerNodeUtil.ByFirstKmerStartPosition);
		
		int weightIn = totalWeight(in);
		assertEquals(2*4+2*5+2*4+2*5, weightIn);
		List<KmerPathNode> list = Lists.newArrayList(new PathSimplificationIterator(in.iterator(), 64, 64));
		assertEquals(1, list.size());
		assertEquals(2, list.get(0).length());
		assertEquals(4, list.get(0).width());
		assertEquals(1, list.get(0).startPosition(0));
		assertEquals(4, list.get(0).endPosition(0));
		assertEquals(4+5, list.get(0).weight());
		assertEquals(weightIn, totalWeight(list));
	}
}
