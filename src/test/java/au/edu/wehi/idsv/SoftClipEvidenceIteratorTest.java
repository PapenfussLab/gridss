package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import htsjdk.samtools.SAMRecord;

import java.util.List;

import org.junit.Before;
import org.junit.Test;

import com.google.common.collect.Lists;


public class SoftClipEvidenceIteratorTest extends TestHelper {
	private List<SAMRecord> sv;
	private List<SoftClipEvidence> out;
	@Before
	public void setup() {
		sv = Lists.newArrayList();
		out = Lists.newArrayList();
	}
	public void go() {
		sv = sorted(sv);
		out = Lists.newArrayList(new SoftClipEvidenceIterator(getContext(), SES(), sv.iterator()));
		// check output is in order
		//for (int i = 0; i < out.size() - 1; i++) {
		//	BreakendSummary l0 = out.get(i).getBreakendSummary();
		//	BreakendSummary l1 = out.get(i).getBreakendSummary();
		//	assertTrue(l0.referenceIndex < l1.referenceIndex || (l0.referenceIndex == l1.referenceIndex && l0.start <= l1.start));
		//}
	}
	@Test
	public void should_return_sc() {
		SAMRecord r = Read(0, 1, "5S10M5S");
		sv.add(r);
		go();
		// forward and backward
		assertEquals(2, out.size());
		assertTrue(out.get(0) instanceof SoftClipEvidence);
		assertTrue(out.get(1) instanceof SoftClipEvidence);
	}
	@Test
	public void should_ignore_filtered_softclips() {
		sv.add(withMapq(0, withReadName("ReadName", Read(0, 1, "5S10M5S")))[0]);
		go();
		assertEquals(0, out.size());
	}
}
