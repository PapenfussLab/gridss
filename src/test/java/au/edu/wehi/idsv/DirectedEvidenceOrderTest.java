package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.stream.Collectors;

import org.junit.Test;

import com.google.common.collect.Lists;

public class DirectedEvidenceOrderTest extends TestHelper {
	@Test
	public void ByNatural_should_be_consistent_with_BreakendSummary_ByStartEnd() {
		MockDirectedEvidence de2 = new MockDirectedEvidence(new BreakendSummary(0, FWD, 2, 1, 3));
		MockDirectedEvidence de3 = new MockDirectedEvidence(new BreakendSummary(0, FWD, 3, 1, 3));
		ArrayList<MockDirectedEvidence> list = Lists.newArrayList(de3, de2);
		assertEquals(BreakendSummary.ByStartEnd.isOrdered(list.stream().map(de -> de.getBreakendSummary()).collect(Collectors.toList())),
				DirectedEvidenceOrder.ByNatural.isOrdered(list));
	}
}
