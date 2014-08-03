package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import htsjdk.samtools.SAMRecord;

import java.util.Collections;
import java.util.List;

import org.junit.Before;
import org.junit.Test;

import com.google.common.collect.Lists;


public class RealignedSoftClipEvidenceIteratorTest extends TestHelper {
	private List<SAMRecord> sv;
	private List<SAMRecord> realigned;	
	private List<SoftClipEvidence> out;
	@Before
	public void setup() {
		sv = Lists.newArrayList();
		realigned = Lists.newArrayList();
		out = Lists.newArrayList();
	}
	public void go() {
		sv = sorted(sv);
		List<SoftClipEvidence> sclist = Lists.newArrayList(new SoftClipEvidenceIterator(getContext(), SES(), sv.iterator()));
		Collections.sort(sclist, DirectedEvidenceOrder.ByNatural);
		out = Lists.newArrayList(new RealignedSoftClipEvidenceIterator(sclist.iterator(), realigned.iterator()));
	}
	@Test
	public void should_allow_realign_in_order_at_same_position() {
		SAMRecord r = withReadName("ReadName", Read(0, 1, "5S10M5S"))[0];
		sv.add(r);
		SAMRecord f = withReadName("0#10#fReadName", Read(0, 1, "5M"))[0];
		SAMRecord b = withReadName("0#1#bReadName", Read(0, 1, "5M"))[0];
		realigned.add(b);
		realigned.add(f);
		go();
		assertEquals(2, out.size());
		assertTrue(out.get(0).getBreakendSummary() instanceof BreakpointSummary);
		assertTrue(out.get(1).getBreakendSummary() instanceof BreakpointSummary);
	}
	@Test
	public void should_require_realign_in_call_position_order() {
		SAMRecord r = withReadName("ReadName", Read(0, 1, "5S10M5S"))[0];
		sv.add(r);
		realigned.add(withReadName("0#1#bReadName", Read(0, 1, "5M"))[0]);
		realigned.add(withReadName("0#10#fReadName", Read(0, 1, "5M"))[0]);
		go();
		assertEquals(2, out.size());
		assertTrue(out.get(0) instanceof RealignedSoftClipEvidence); // backward
		assertTrue(out.get(1) instanceof RealignedSoftClipEvidence); // forward
	}
}
