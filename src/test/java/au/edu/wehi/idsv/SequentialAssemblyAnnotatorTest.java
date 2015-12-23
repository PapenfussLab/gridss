package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;

import java.util.List;

import org.junit.Test;

import com.google.common.collect.Lists;


public class SequentialAssemblyAnnotatorTest extends TestHelper {
	@Test
	public void should_hydrate_indel_parent() {
		DirectedEvidence e0 = SCE(FWD, Read(0, 1, "1M12S"));
		SAMRecordAssemblyEvidence a1 = AssemblyFactory.createAnchoredBreakpoint(getContext(), AES(), Lists.transform(Lists.newArrayList(e0), EID), 0, 1, 1, 0, 2, 1, B("ATA"), B("ATA"));
		List<SAMRecordAssemblyEvidence> result = Lists.newArrayList(new SequentialAssemblyAnnotator(getContext().getLinear(),
				Lists.newArrayList(a1).iterator(),
				Lists.newArrayList(e0).iterator(),
				1));
		assertEquals(1, result.size());
		assertEquals(1, result.get(0).getAssemblySupportCount());
	}
	@Test
	public void should_hydrate_all_assemblies() {
		DirectedEvidence e0 = SCE(FWD, Read(0, 1, "1M12S"));
		DirectedEvidence e1 = SCE(FWD, Read(0, 10, "1M10S"));
		DirectedEvidence e2 = SCE(FWD, Read(0, 10, "1M11S"));
		DirectedEvidence e3 = SCE(FWD, Read(0, 10, "1M12S"));
		DirectedEvidence e4 = SCE(FWD, Read(0, 100, "1M12S"));
		SAMRecordAssemblyEvidence a1 = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), FWD, Lists.transform(Lists.newArrayList(e2, e3), EID), 0, 10, 1, B("TT"), B("TT"));
		SAMRecordAssemblyEvidence a2 = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), FWD, Lists.transform(Lists.newArrayList(e1, e2), EID), 0, 11, 1, B("TT"), B("TT"));
		List<SAMRecordAssemblyEvidence> result = Lists.newArrayList(new SequentialAssemblyAnnotator(getContext().getLinear(),
				Lists.newArrayList(a1, a2).iterator(),
				Lists.newArrayList(e0, e1, e2, e3, e4).iterator(),
				10));
		assertEquals(2, result.size());
		assertEquals(2, result.get(0).getAssemblySupportCount());
		assertEquals(2, result.get(1).getAssemblySupportCount());
	}
}
