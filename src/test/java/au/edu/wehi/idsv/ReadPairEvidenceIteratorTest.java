package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import htsjdk.samtools.SAMRecord;

import java.util.List;

import org.junit.Before;
import org.junit.Test;

import com.google.common.collect.Lists;


public class ReadPairEvidenceIteratorTest extends TestHelper {
	private List<SAMRecord> sv;
	private List<SAMRecord> mate;
	private List<NonReferenceReadPair> out;
	@Before
	public void setup() {
		sv = Lists.newArrayList();
		mate = Lists.newArrayList();
		out = Lists.newArrayList();
	}
	public void go() {
		sv = sorted(sv);
		mate = mateSorted(mate);
		ReadPairEvidenceIterator it = new ReadPairEvidenceIterator(
				getContext(),
				SES(),
				sv.iterator(),
				mate.iterator());
		out = Lists.newArrayList(it);
		// check output is in order
		for (int i = 0; i < out.size() - 1; i++) {
			BreakendSummary l0 = out.get(i).getBreakendSummary();
			BreakendSummary l1 = out.get(i).getBreakendSummary();
			assertTrue(l0.referenceIndex < l1.referenceIndex || (l0.referenceIndex == l1.referenceIndex && l0.start <= l1.start));
		}
	}
	@Test
	public void should_pair_oea_with_mate() {
		sv.add(OEA(0, 1, "100M", true)[0]);
		mate.add(OEA(0, 1, "100M", true)[1]);
		go();
		assertEquals(1, out.size());
		assertTrue(out.get(0) instanceof NonReferenceReadPair);
	}
	@Test
	public void should_pair_dp_with_mate() {
		sv.add(DP(0, 1, "100M", true, 1, 1, "100M", true)[0]);
		mate.add(DP(0, 1, "100M", true, 1, 1, "100M", true)[1]);
		go();
		assertEquals(1, out.size());
		assertTrue(out.get(0) instanceof NonReferenceReadPair);
	}
	@Test
	public void should_ignore_non_sv_reads() {
		sv.add(RP(0, 1, 2, 1)[0]);
		sv.add(RP(0, 1, 2, 1)[1]);
		go();
		assertEquals(0, out.size());
	}
	@Test
	public void should_expect_mates_ordered_by_mate_alignment_start() {
		sv.add(withReadName("DP", DP(0, 2, "100M", true, 1, 1, "100M", true))[0]);
		mate.add(withReadName("DP", DP(0, 2, "100M", true, 1, 1, "100M", true))[1]);
		sv.add(withReadName("OEA", OEA(0, 1, "100M", true))[0]);
		mate.add(withReadName("OEA", OEA(0, 1, "100M", true))[1]);
		go();
		assertEquals(2, out.size());
		assertTrue(out.get(0) instanceof NonReferenceReadPair);
		assertTrue(out.get(1) instanceof NonReferenceReadPair);
	}
	@Test
	public void should_filter_if_local_mapq_too_low() {
		sv.add(withMapq(1, OEA(0, 1, "100M", true))[0]);
		mate.add(withMapq(44, OEA(1, 1, "100M", true))[1]);
		go();
		assertEquals(0, out.size());
	}
	@Test
	public void should_not_filter_on_remote_mapq() {
		sv.add(withMapq(44, DP(0, 1, "5M", true, 0, 1, "5M", true))[0]);
		mate.add(withMapq(1, DP(0, 1, "5M", true, 0, 1, "5M", true))[1]);
		go();
		assertEquals(1, out.size());
		assertTrue(out.get(0) instanceof NonReferenceReadPair);
	}
	@Test
	public void should_filter_self_overlap() {
		// base 5 common to both
		sv.add(DP(0, 1, "5M", true, 0, 5, "5M", false)[0]);
		mate.add(DP(0, 1, "5M", true, 0, 5, "5M", false)[1]);
		// adjacent but not overlapping
		sv.add(DP(1, 1, "5M", true, 1, 6, "5M", false)[0]);
		mate.add(DP(1, 1, "5M", true, 1, 6, "5M", false)[1]);
		go();
		assertEquals(1, out.size());
	}
	@Test
	public void should_not_filter_self_overlap_in_unexpected_orientation() {
		// overlapping unexpected orientation
		sv.add(DP(0, 1, "5M", true, 0, 5, "5M", true)[0]);
		mate.add(DP(0, 1, "5M", true, 0, 5, "5M", true)[1]);
		go();
		assertEquals(1, out.size());
	}
}
