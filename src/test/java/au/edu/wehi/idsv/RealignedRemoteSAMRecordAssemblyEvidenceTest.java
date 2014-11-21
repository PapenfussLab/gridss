package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;

import org.junit.Test;


public class RealignedRemoteSAMRecordAssemblyEvidenceTest extends TestHelper {
	@Test
	public void should_swap_breakend_position() {
		BreakpointSummary bs = new BreakpointSummary(0, FWD, 1, 1, 1, BWD, 3, 3);
		SAMRecordAssemblyEvidence e = A(bs); 
		RealignedRemoteSAMRecordAssemblyEvidence r = new RealignedRemoteSAMRecordAssemblyEvidence(getContext(), AES(), e.getSAMRecord(), e.getRemoteSAMRecord());
		assertEquals(r.getBreakendSummary(), bs.remoteBreakpoint());
	}
}
