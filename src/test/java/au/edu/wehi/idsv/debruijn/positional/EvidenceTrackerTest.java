package au.edu.wehi.idsv.debruijn.positional;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;

import com.google.common.collect.ImmutableList;


public class EvidenceTrackerTest extends TestHelper {
	@Test
	public void shouldTrackWhenIterated() {
		int k = 4;
		List<KmerSupportNode> list = new ArrayList<KmerSupportNode>();
		list.add(KmerEvidence.create(k, SCE(FWD, Read(0, 1, "4M1S")), false).node(0));
		EvidenceTracker tracker = new EvidenceTracker();
		Set<KmerEvidence> result = tracker.untrack(ImmutableList.of(new KmerPathSubnode(KPN(k, "AAAA", 1, 1, true))));
		assertEquals(0, result.size());
		tracker.track(list.get(0));
		result = tracker.untrack(ImmutableList.of(new KmerPathSubnode(KPN(k, "AAAA", 1, 1, true))));
		assertEquals(1, result.size());
		assertEquals(list.get(0).evidence(), result.iterator().next());
	}
	@Test
	public void shouldUntrackAllEvidenceNodes() {
		int k = 4;
		KmerEvidence e = KmerEvidence.create(k, SCE(FWD, Read(0, 1, "4M1S")), false);
		List<KmerSupportNode> list = new ArrayList<KmerSupportNode>();
		list.add(e.node(0));
		list.add(e.node(1));
		EvidenceTracker tracker = new EvidenceTracker();
		tracker.track(list.get(0));
		tracker.track(list.get(1));
		Set<KmerEvidence> result = tracker.untrack(ImmutableList.of(new KmerPathSubnode(KPN(k, "AAAA", 1, 1, true))));
		assertEquals(1, result.size());
		result = tracker.untrack(ImmutableList.of(new KmerPathSubnode(KPN(k, "AAAA", 1, 1, true))));
		assertEquals(0, result.size());
		result = tracker.untrack(ImmutableList.of(new KmerPathSubnode(KPN(k, "AAAA", 2, 2, true))));
		assertEquals(0, result.size());
	}
	@Test
	public void shouldTrackEvidenceID() {
		int k = 4;
		KmerEvidence e = KmerEvidence.create(k, SCE(FWD, Read(0, 1, "4M1S")), false);
		List<KmerSupportNode> list = new ArrayList<KmerSupportNode>();
		list.add(e.node(0));
		EvidenceTracker tracker = new EvidenceTracker();
		assertFalse(tracker.isTracked(e.evidenceId()));
		tracker.track(list.get(0));
		assertTrue(tracker.isTracked(e.evidenceId()));
		tracker.remove(e);
		assertFalse(tracker.isTracked(e.evidenceId()));
	}
}
