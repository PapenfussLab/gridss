package au.edu.wehi.idsv.debruijn.positional;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

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
		EvidenceTracker tracker = new EvidenceTracker(list.iterator());
		Set<KmerEvidence> result = tracker.untrack(ImmutableList.of(new KmerPathSubnode(KPN(k, "AAAA", 1, 1, true))));
		assertEquals(0, result.size());
		assertEquals(list.get(0), tracker.next());
		assertFalse(tracker.hasNext());
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
		EvidenceTracker tracker = new EvidenceTracker(list.iterator());
		tracker.next();
		tracker.next();
		Set<KmerEvidence> result = tracker.untrack(ImmutableList.of(new KmerPathSubnode(KPN(k, "AAAA", 1, 1, true))));
		assertEquals(1, result.size());
		result = tracker.untrack(ImmutableList.of(new KmerPathSubnode(KPN(k, "AAAA", 1, 1, true))));
		assertEquals(0, result.size());
		result = tracker.untrack(ImmutableList.of(new KmerPathSubnode(KPN(k, "AAAA", 2, 2, true))));
		assertEquals(0, result.size());
	}
}
