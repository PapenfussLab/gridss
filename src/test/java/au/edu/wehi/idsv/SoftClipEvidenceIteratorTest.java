package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import htsjdk.samtools.SAMRecord;

import java.util.List;
import java.util.stream.Collectors;

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
		out = Lists.newArrayList(new SoftClipEvidenceIterator(SES(), sv.iterator()));
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
	@Test
	public void should_ignore_unmapped_softclips() {
		SAMRecord r = Read(0, 1, "50S50M");
		r.setReadUnmappedFlag(true);
		sv.add(r);
		go();
		assertEquals(0, out.size());
	}
	@Test
	public void should_treat_indels_as_realigned_soft_clips() {
		SAMRecord r = Read(1, 1, "2M1D2M");
		sv.add(r);
		go();
		assertEquals(2, out.size());
		assertEquals(new BreakpointSummary(1, FWD, 2, 2, 1, BWD, 4, 4), out.get(0).getBreakendSummary());
		assertEquals(new BreakpointSummary(1, BWD, 4, 4, 1, FWD, 2, 2), out.get(1).getBreakendSummary());
	}
	@Test
	public void should_extract_all_indels() {
		SAMRecord r = Read(1, 1, "2M1D2M1D1M");
		sv.add(r);
		go();
		assertEquals(4, out.size());
	}
	@Test
	public void should_have_unique_names() {
		SAMRecord r = Read(1, 1, "5S5M5I5M5D5M5S");
		sv.add(r);
		go();
		assertEquals(out.size(), out.stream().map(sce -> sce.getEvidenceID()).collect(Collectors.toSet()).size());
	}
}
