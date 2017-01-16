package au.edu.wehi.idsv;

import java.io.IOException;
import java.util.Random;

import org.junit.Ignore;
import org.junit.Test;

import au.edu.wehi.idsv.GreedyAllocationCache.EventAlignmentScoreNode;
import au.edu.wehi.idsv.GreedyAllocationCache.Hash96bit;
import htsjdk.samtools.util.Log;

public class GreedyAllocationCacheTest {
	private static final Log log = Log.getInstance(GreedyAllocationCacheTest.class);
	public class TestCache extends GreedyAllocationCache {
		public GreedyAllocationCacheLookup<EventAlignmentScoreNode> cache = createLookup("test", EventAlignmentScoreNode.class, 100000000);
		@Override
		public void close() throws IOException {
			cache.close();
		}
	}
	@Test
	@Ignore("Too much memory for desktop testing")
	public void should_use_off_heap_memory() throws IOException {
		Random rnd = new Random();
		try (TestCache cache = new TestCache()) {
			for (int i = 0; i < Integer.MAX_VALUE; i++) {
				cache.cache.put(new Hash96bit(rnd.nextInt(), rnd.nextInt()),
						new EventAlignmentScoreNode(new Hash96bit(rnd.nextInt(), rnd.nextInt()),
								rnd.nextInt(10000),
								new Hash96bit(rnd.nextInt(), rnd.nextInt())));
				if (i % 100000 == 0) {
					log.warn(String.format("Loaded %,d records. Current java heap memory usage is %,d MiB",
							i,
							(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) >> 20));
				}
			}
		}
	}
}
