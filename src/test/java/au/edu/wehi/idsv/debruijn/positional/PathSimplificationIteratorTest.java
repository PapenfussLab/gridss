package au.edu.wehi.idsv.debruijn.positional;

import static org.junit.Assert.*;

import java.util.List;

import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;


public class PathSimplificationIteratorTest extends TestHelper {
	/**
	 * Create a data set that tiles
	 */
	List<KmerNode> tiled(int width, int length, int n) {
		int k = 25;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				// Need to tile in multiple steps
				new KmerAggregateNode(kmer, weight, start, end, reference)
				ImmutableKmerNode node = new ImmutableKmerNode(kmer, start, end, weight, isReference);
				Evidence e = new Evidence(String.format("i=%d,j=%d", i, j), i * width, error, k, bases, new byte[bases.length], false, false)
			}
		}
	}
	KmerPathNode KPN(int start, int end, int length) {
		KmerPathNode pn = new KmerPathNode(0, 1, start, end, false);
		for (int i = 1; i < length; i++) {
			pn.append(new KmerAggregateNode(0, 1, start + i, end + i, false));
		}
		return pn;
	}
	@Test
	public void test_worst_case_tiling() {
		for (int maxSupportWidth = 1; maxSupportWidth < 16; maxSupportWidth++) {
			for (int maxPathLength = 1; maxPathLength < 16; maxPathLength++) {
				for (int supportWidth = 1; supportWidth <= maxSupportWidth; supportWidth++) {
					for (int pathLength = 1; pathLength <= maxPathLength; pathLength++) {
						// tile evidence such that the entire input could be compressed to a single
						// node if the limits did not exist
						List<KmerPathNode> input = 
						new KmerPathNode(kmer, weight, start, end, reference)
						fail();
					}
				}
			}
		}
	}
}
