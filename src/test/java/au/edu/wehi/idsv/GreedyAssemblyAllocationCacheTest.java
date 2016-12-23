package au.edu.wehi.idsv;

import static org.junit.Assert.*;

import org.junit.Test;

import com.google.common.collect.ImmutableList;

public class GreedyAssemblyAllocationCacheTest extends TestHelper {
	@Test
	public void should_only_consider_local_read_mapping() {
		GreedyAssemblyAllocationCache cache = new GreedyAssemblyAllocationCache();
		SoftClipEvidence r1 = SCE(FWD, Read(0, 1, "1M1S"));
		SoftClipEvidence r2 = SCE(FWD, Read(0, 2, "1M1S"));
		cache.addBreakendAssemblyAllocation(1, r1);
		cache.addBreakendAssemblyAllocation(1, r2);
		assertTrue(cache.isBestBreakendAssemblyAllocation(r1));
		assertTrue(cache.isBestBreakendAssemblyAllocation(r2));
	}
	@Test
	public void should_merge_allocations() {
		SoftClipEvidence r1a = SCE(FWD, withReadName("r1", Read(0, 1, "1M1S"))[0]);
		SoftClipEvidence r2a = SCE(FWD, withReadName("r2", Read(0, 2, "1M1S"))[0]);
		SoftClipEvidence r1b = SCE(FWD, withReadName("r1", Read(1, 1, "1M1S"))[0]);
		SoftClipEvidence r2b = SCE(FWD, withReadName("r2", Read(1, 2, "1M1S"))[0]);
		GreedyAssemblyAllocationCache cachea = new GreedyAssemblyAllocationCache();
		GreedyAssemblyAllocationCache cacheb = new GreedyAssemblyAllocationCache();
		cachea.addBreakendAssemblyAllocation(1, r1a);
		cachea.addBreakendAssemblyAllocation(2, r2a);
		cacheb.addBreakendAssemblyAllocation(2, r1b);
		cacheb.addBreakendAssemblyAllocation(1, r2b);
		GreedyAssemblyAllocationCache merged = GreedyAssemblyAllocationCache.merge(ImmutableList.of(cachea, cacheb));
		assertFalse(merged.isBestBreakendAssemblyAllocation(r1a));
		assertTrue(merged.isBestBreakendAssemblyAllocation(r1b));
		assertTrue(merged.isBestBreakendAssemblyAllocation(r2a));
		assertFalse(merged.isBestBreakendAssemblyAllocation(r2b));
	}
}
