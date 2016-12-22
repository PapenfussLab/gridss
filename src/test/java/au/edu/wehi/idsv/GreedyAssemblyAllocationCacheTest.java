package au.edu.wehi.idsv;

import static org.junit.Assert.assertTrue;

import org.junit.Test;

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
}
