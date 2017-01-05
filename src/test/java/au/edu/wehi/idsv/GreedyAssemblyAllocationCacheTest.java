package au.edu.wehi.idsv;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.IOException;

import org.junit.Test;

public class GreedyAssemblyAllocationCacheTest extends TestHelper {
	@Test
	public void should_only_consider_local_read_mapping() throws IOException {
		try (GreedyAssemblyAllocationCache cache = new GreedyAssemblyAllocationCache(1000)) {
			SoftClipEvidence r1 = SCE(FWD, Read(0, 1, "1M1S"));
			SoftClipEvidence r2 = SCE(FWD, Read(0, 2, "1M1S"));
			cache.addBreakendAssemblyAllocation(1, r1);
			cache.addBreakendAssemblyAllocation(1, r2);
			assertTrue(cache.isBestBreakendAssemblyAllocation(r1));
			assertTrue(cache.isBestBreakendAssemblyAllocation(r2));
		}
	}
	@Test
	public void should_merge_allocations() throws IOException {
		SoftClipEvidence r1a = SCE(FWD, withReadName("r1", Read(0, 1, "1M1S"))[0]);
		SoftClipEvidence r2a = SCE(FWD, withReadName("r2", Read(0, 2, "1M1S"))[0]);
		SoftClipEvidence r1b = SCE(FWD, withReadName("r1", Read(1, 1, "1M1S"))[0]);
		SoftClipEvidence r2b = SCE(FWD, withReadName("r2", Read(1, 2, "1M1S"))[0]);
		try (GreedyAssemblyAllocationCache cachea = new GreedyAssemblyAllocationCache(1000)) {
			cachea.addBreakendAssemblyAllocation(1, r1a);
			cachea.addBreakendAssemblyAllocation(2, r2a);
			cachea.addBreakendAssemblyAllocation(2, r1b);
			cachea.addBreakendAssemblyAllocation(1, r2b);
			GreedyAssemblyAllocationCache merged = cachea;
			assertFalse(merged.isBestBreakendAssemblyAllocation(r1a));
			assertTrue(merged.isBestBreakendAssemblyAllocation(r1b));
			assertTrue(merged.isBestBreakendAssemblyAllocation(r2a));
			assertFalse(merged.isBestBreakendAssemblyAllocation(r2b));
		}
	}
}
