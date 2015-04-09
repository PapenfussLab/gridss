package au.edu.wehi.idsv.visualisation;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;


public class SubgraphAlgorithmMetricsTest extends TestHelper {
	@Test
	public void should_use_GFF3_compliant_attributes() {
		SubgraphAlgorithmMetrics m = new SubgraphAlgorithmMetrics(getContext(), 0);
		assertFalse(m.toBed().replace("; ", "_").contains(";"));
	}
	@Test
	public void should_write_0_time_for_skipped_steps() {
		SubgraphAlgorithmMetrics m = new SubgraphAlgorithmMetrics(getContext(), 0);
		m.assemblyStarted();
		m.assemblyComplete();
		String s = m.toBed();
		String times = s.substring(s.indexOf("Times") + 5, s.indexOf(";", s.indexOf("Times")) - 1);
		for (String ss : times.split(",")) {
			assertTrue(new Integer(ss.trim()).intValue() >= 0);
		}
	}
	
}
