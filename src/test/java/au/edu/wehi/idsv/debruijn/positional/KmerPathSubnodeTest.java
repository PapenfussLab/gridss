package au.edu.wehi.idsv.debruijn.positional;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.junit.Test;

import com.google.common.collect.Lists;
import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;

import au.edu.wehi.idsv.TestHelper;


public class KmerPathSubnodeTest extends TestHelper {
	@Test
	public void givenNext_should_restrict_range() {
		KmerPathNode pn = KPN(4, "AAAA", 1, 10, true);
		assertSame(new KmerPathSubnode(pn, 4, 6), new KmerPathSubnode(pn).givenNext(new KmerPathSubnode(KPN(4, "AAAT", 5, 7, true))));
		assertSame(new KmerPathSubnode(pn, 1, 10), new KmerPathSubnode(pn).givenNext(new KmerPathSubnode(KPN(4, "AAAT", -100, 100, true))));
		assertSame(new KmerPathSubnode(pn, 1, 10), new KmerPathSubnode(pn).givenNext(new KmerPathSubnode(KPN(4, "AAAT", 1, 100, true))));
		assertSame(new KmerPathSubnode(pn, 1, 10), new KmerPathSubnode(pn).givenNext(new KmerPathSubnode(KPN(4, "AAAT", 2, 100, true))));
		assertSame(new KmerPathSubnode(pn, 2, 10), new KmerPathSubnode(pn).givenNext(new KmerPathSubnode(KPN(4, "AAAT", 3, 100, true))));
		assertSame(new KmerPathSubnode(pn, 2, 10), new KmerPathSubnode(pn).givenNext(new KmerPathSubnode(KPN(4, "AAAT", 3, 12, true))));
		assertSame(new KmerPathSubnode(pn, 2, 10), new KmerPathSubnode(pn).givenNext(new KmerPathSubnode(KPN(4, "AAAT", 3, 11, true))));
		assertSame(new KmerPathSubnode(pn, 2, 9), new KmerPathSubnode(pn).givenNext(new KmerPathSubnode(KPN(4, "AAAT", 3, 10, true))));
	}
	@Test
	public void givenPrev_should_restrict_range() {
		KmerPathNode pn = KPN(4, "AAAA", 1, 10, true);
		assertSame(new KmerPathSubnode(pn, 6, 8), new KmerPathSubnode(pn).givenPrev(new KmerPathSubnode(KPN(4, "TAAA", 5, 7, true))));
		assertSame(new KmerPathSubnode(pn, 1, 10), new KmerPathSubnode(pn).givenPrev(new KmerPathSubnode(KPN(4, "TAAA", -100, 100, true))));
		assertSame(new KmerPathSubnode(pn, 1, 10), new KmerPathSubnode(pn).givenPrev(new KmerPathSubnode(KPN(4, "TAAA", -1, 100, true))));
		assertSame(new KmerPathSubnode(pn, 1, 10), new KmerPathSubnode(pn).givenPrev(new KmerPathSubnode(KPN(4, "TAAA", 0, 100, true))));
		assertSame(new KmerPathSubnode(pn, 2, 10), new KmerPathSubnode(pn).givenPrev(new KmerPathSubnode(KPN(4, "TAAA", 1, 100, true))));
		assertSame(new KmerPathSubnode(pn, 4, 10), new KmerPathSubnode(pn).givenPrev(new KmerPathSubnode(KPN(4, "TAAA", 3, 10, true))));
		assertSame(new KmerPathSubnode(pn, 4, 10), new KmerPathSubnode(pn).givenPrev(new KmerPathSubnode(KPN(4, "TAAA", 3, 9, true))));
		assertSame(new KmerPathSubnode(pn, 4, 9), new KmerPathSubnode(pn).givenPrev(new KmerPathSubnode(KPN(4, "TAAA", 3, 8, true))));
	}
	@Test
	public void next_should_restrict_range() {
		KmerPathNode pn = KPN(4, "AAAA", 1, 10, true);
		KmerPathNode.addEdge(pn, KPN(4, "AAAA", -100, 100, true));
		KmerPathNode.addEdge(pn, KPN(4, "AAAT", 3, 11, true));
		KmerPathNode.addEdge(pn, KPN(4, "AAAG", 2, 10, true));
		KmerPathNode.addEdge(pn, KPN(4, "AAAC", 4, 5, true));
		KmerPathNode.addEdge(pn, KPN(4, "AAAC", 6, 6, true));
		assertEquals(5, new KmerPathSubnode(pn).next().size());
		assertEquals(2, new KmerPathSubnode(pn, 1, 1).next().size());
		assertSame(new KmerPathSubnode(KPN(4, "AAAA", -100, 100, true), 2, 2), new KmerPathSubnode(pn, 1, 1).next().get(0));
	}
	@Test
	public void prev_should_restrict_range() {
		KmerPathNode pn = KPN(4, "AAAA", 1, 10, true);
		KmerPathNode.addEdge(KPN(4, "AAAA", -100, 100, true), pn);
		KmerPathNode.addEdge(KPN(4, "TAAA", 1, 9, true), pn);
		KmerPathNode.addEdge(KPN(4, "GAAA", 0, 8, true), pn);
		KmerPathNode.addEdge(KPN(4, "CAAA", 2, 3, true), pn);
		KmerPathNode.addEdge(KPN(4, "CAAA", 4, 4, true), pn);
		assertEquals(5, new KmerPathSubnode(pn).prev().size());
		assertEquals(2, new KmerPathSubnode(pn, 1, 1).prev().size());
		assertSame(new KmerPathSubnode(KPN(4, "AAAA", -100, 100, true), 0, 0), new KmerPathSubnode(pn, 1, 1).prev().get(0));
	}
	@Test
	public void nextPathRangesOfDegree() {
		//  1234567890
		//  2344433332 degree
		//------------
		//  ---------
		//   ---------
		//    --
		//      -
		KmerPathNode pn = KPN(4, "AAAA", 1, 10, true);
		KmerPathNode.addEdge(pn, KPN(4, "AAAA", -100, 100, true));
		KmerPathNode.addEdge(pn, KPN(4, "AAAT", 3, 11, true));
		KmerPathNode.addEdge(pn, KPN(4, "AAAG", 2, 10, true));
		KmerPathNode.addEdge(pn, KPN(4, "AAAC", 4, 5, true));
		KmerPathNode.addEdge(pn, KPN(4, "AAAC", 6, 6, true));
		
		assertEquals(TreeRangeSet.create(), new KmerPathSubnode(pn).nextPathRangesOfDegree(0));
		
		assertEquals(TreeRangeSet.create(), new KmerPathSubnode(pn).nextPathRangesOfDegree(1));
		
		RangeSet<Integer> expected2 = TreeRangeSet.create();
		expected2.add(Range.closed(1, 1));
		expected2.add(Range.closed(10, 10));
		assertEquals(expected2, new KmerPathSubnode(pn).nextPathRangesOfDegree(2));
		
		RangeSet<Integer> expected3 = TreeRangeSet.create();
		expected3.add(Range.closed(2, 2));
		expected3.add(Range.closed(6, 9));
		assertEquals(expected3, new KmerPathSubnode(pn).nextPathRangesOfDegree(3));
		
		RangeSet<Integer> expected4 = TreeRangeSet.create();
		expected4.add(Range.closed(3, 4));
		expected4.add(Range.closed(5, 5));
		assertEquals(expected4, new KmerPathSubnode(pn).nextPathRangesOfDegree(4));
	
		assertEquals(TreeRangeSet.create(), new KmerPathSubnode(pn).nextPathRangesOfDegree(5));
	}
	@Test
	public void prevPathRangesOfDegree() {
		KmerPathNode pn = KPN(4, "AAAA", -100, 100, true);
		KmerPathNode.addEdge(KPN(4, "AAAA", -100, 9, true), pn);
		KmerPathNode.addEdge(KPN(4, "AAAA", 9, 9, true), pn);
		KmerPathNode.addEdge(KPN(4, "AAAA", 15, 15, true), pn);
		KmerPathNode.addEdge(KPN(4, "AAAA", 15, 16, true), pn);
		KmerPathNode.addEdge(KPN(4, "AAAA", 15, 17, true), pn);
		KmerPathNode.addEdge(KPN(4, "AAAAA", 15, 15, true), pn);
		KmerPathNode.addEdge(KPN(4, "AAAA", 50, 100, true), pn);
		//         1         2
		//12345678901234567890
		//----------     -
		//----------     --
		//               ---
		//                -
		KmerPathSubnode sn = new KmerPathSubnode(pn, 10, 20);
		
		RangeSet<Integer> expected0 = TreeRangeSet.create();
		expected0.add(Range.closed(11, 15));
		expected0.add(Range.closed(19, 20));
		assertEquals(expected0, sn.prevPathRangesOfDegree(0));
		
		RangeSet<Integer> expected1 = TreeRangeSet.create();
		expected1.add(Range.closed(18, 18));
		assertEquals(expected1, sn.prevPathRangesOfDegree(1));
		
		RangeSet<Integer> expected2 = TreeRangeSet.create();
		expected2.add(Range.closed(10, 10));
		assertEquals(expected2, sn.prevPathRangesOfDegree(2));
		
		RangeSet<Integer> expected3 = TreeRangeSet.create();
		expected3.add(Range.closed(16, 16));
		expected3.add(Range.closed(17, 17));
		assertEquals(expected3, sn.prevPathRangesOfDegree(3));
		
		RangeSet<Integer> expected4 = TreeRangeSet.create();
		assertEquals(expected4, sn.prevPathRangesOfDegree(4));
	}
	@Test
	public void width_should_match_firstKmer_width() {
		KmerPathNode pn = KPN(4, "AAAA", 1, 10, true);
		assertEquals(2, new KmerPathSubnode(pn, 2, 3).width());
	}
	@Test
	public void next_nextPathRangesOfDegree_regression_test1() {
		int k = 4;
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(k, "TAAAAA", 1, 10, false));
		input.add(KPN(k, "AAATC", 5, 5, false));
		input.add(KPN(k, "AAAGC", 8, 9, false));
		KmerPathNode.addEdge(input.get(0), input.get(1));
		KmerPathNode.addEdge(input.get(0), input.get(2));
		KmerPathSubnode sn = new KmerPathSubnode(input.get(0));
		
		assertEquals(2, sn.next().size());
		
		RangeSet<Integer> expected = TreeRangeSet.create();
		expected.add(Range.closed(1, 1));
		expected.add(Range.closed(3, 4));
		expected.add(Range.closed(7, 10));
		assertEquals(expected, sn.nextPathRangesOfDegree(0));
	}
	@Test
	public void should_traverse_simple_sc() { 
		List<KmerPathNode> paths = Lists.newArrayList(asKPN(4, 100, SCE(FWD, withSequence("ACGTGGTCGACC", Read(0, 5, "6M6S")))));
		assertEquals(2, paths.size());
		assertEquals(1, new KmerPathSubnode(paths.get(0)).next().size());
		assertEquals(1, new KmerPathSubnode(paths.get(1)).prev().size());
	}
	@Test
	public void subnodesOfDegree() {
		int k = 4;
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(k, "TAAA", 0, 9, false));
		input.add(KPN(k, "CAAA", 4, 6, false));
		input.add(KPN(k, "AAAA", 1, 10, false));
		input.add(KPN(k, "AAAC", 5, 9, false));
		input.add(KPN(k, "AAAT", 10, 10, false));
		input.add(KPN(k, "AAAG", 4, 5, false));
		KmerPathNode.addEdge(input.get(0), input.get(2));
		KmerPathNode.addEdge(input.get(1), input.get(2));
		KmerPathNode.addEdge(input.get(2), input.get(3));
		KmerPathNode.addEdge(input.get(2), input.get(4));
		KmerPathNode.addEdge(input.get(2), input.get(5));
		List<KmerPathSubnode> result = new KmerPathSubnode(input.get(2)).subnodesOfDegree(2, 1);
		// 1111222111 prev
		// ----------
		//     456
		// 0123456789
		// 1234567890 (2)
		//    56789
		//         0
		//   45
		// ----------
		// 0012111110 next
		// 1111222111 (prev)
		// 1234567890
		assertEquals(1, result.size());
		assertEquals(5, result.get(0).firstStart());
		assertEquals(7, result.get(0).firstEnd());
		result = new KmerPathSubnode(input.get(2)).subnodesOfDegree(1, 1);
		assertEquals(3, result.size());
		assertEquals(3, result.get(0).firstStart());
		assertEquals(3, result.get(0).firstEnd());
		// {[3,3] [8,8] [9,9]} not  {[3,3] [8,9]} since the adj nodes are different for 8 and 9
		assertEquals(8, result.get(1).firstStart());
		assertEquals(8, result.get(1).firstEnd());
		assertEquals(9, result.get(2).firstStart());
		assertEquals(9, result.get(2).firstEnd());
	}
}
