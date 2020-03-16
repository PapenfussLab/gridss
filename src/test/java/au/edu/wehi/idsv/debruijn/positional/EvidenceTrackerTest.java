package au.edu.wehi.idsv.debruijn.positional;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.NonReferenceReadPair;
import au.edu.wehi.idsv.SingleReadEvidence;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import htsjdk.samtools.SAMRecord;
import org.junit.Test;

import com.google.common.collect.ImmutableList;

import au.edu.wehi.idsv.TestHelper;


public class EvidenceTrackerTest extends TestHelper {
	@Test
	public void shouldTrackWhenIterated() {
		int k = 4;
		List<KmerSupportNode> list = new ArrayList<KmerSupportNode>();
		list.add(KmerEvidence.create(k, SCE(FWD, Read(0, 1, "4M1S"))).node(0));
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
		KmerEvidence e = KmerEvidence.create(k, SCE(FWD, Read(0, 1, "4M1S")));
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
		KmerEvidence e = KmerEvidence.create(k, SCE(FWD, Read(0, 1, "4M1S")));
		List<KmerSupportNode> list = new ArrayList<KmerSupportNode>();
		list.add(e.node(0));
		EvidenceTracker tracker = new EvidenceTracker();
		assertFalse(tracker.isTracked(e.evidence().getEvidenceID()));
		tracker.track(list.get(0));
		assertTrue(tracker.isTracked(e.evidence().getEvidenceID()));
		tracker.remove(Collections.singleton(e));
		assertFalse(tracker.isTracked(e.evidence().getEvidenceID()));
	}
	@Test
	public void should_remove_efficiently_in_degenerate_sequence() {
		final int k = 4;
		final MockSAMEvidenceSource ses = SES();
		Set<KmerEvidence> list = IntStream.range(1, 100000).mapToObj(i ->
			KmerEvidence.create(k, SCE(FWD, ses, withReadName("read" + i, Read(0, 1, "7M1S"))))
		).collect(Collectors.toSet());
		EvidenceTracker tracker = new EvidenceTracker();
		for (int i = 0; i < 8 - k + 1; i++) {
			int finalI = i;
			list.stream().forEach(ke -> tracker.track(ke.node(finalI)));
		}
		tracker.remove(list);
	}
	@Test
	public void remove_should_remove_both_anchor_and_nonreference_KmerEvidence_for_read_pair() {
		int k = 4;
		SAMRecord[] rp = OEA(3, 1, "4M", true);
		rp[0].setReadBases(B("AACC"));
		rp[1].setReadBases(B("ACGT"));
		NonReferenceReadPair nrrp = NonReferenceReadPair.create(rp[0], rp[1], SES(50));
		KmerEvidence e = KmerEvidence.create(k, nrrp);
		KmerEvidence e2 = KmerEvidence.createAnchor(k, nrrp, 0, nrrp.getEvidenceSource().getContext().getReference());
		EvidenceTracker tracker = new EvidenceTracker();
		Stream.of(e, e2)
				.flatMap(ev -> IntStream.range(0, ev.length()).mapToObj(i -> ev.node(i)))
				.forEach(ksn -> tracker.track(ksn));
		assertEquals(2, tracker.getTrackedEvidence().size());
		assertTrue(tracker.isTracked(nrrp.getEvidenceID()));
		tracker.sanityCheck();
		tracker.remove(Collections.singleton(e));
		assertEquals(0, tracker.getTrackedEvidence().size());
		assertFalse(tracker.isTracked(nrrp.getEvidenceID()));
		tracker.sanityCheck();
	}
	@Test
	public void untrack_should_remove_both_anchor_and_nonreference_KmerEvidence_for_read_pair() {
		int k = 4;
		SAMRecord[] rp = OEA(3, 1, "4M", true);
		rp[0].setReadBases(B("AACC"));
		rp[1].setReadBases(B("ACGT"));
		NonReferenceReadPair nrrp = NonReferenceReadPair.create(rp[0], rp[1], SES(50));
		KmerEvidence e = KmerEvidence.create(k, nrrp);
		KmerEvidence e2 = KmerEvidence.createAnchor(k, nrrp, 0, nrrp.getEvidenceSource().getContext().getReference());
		EvidenceTracker tracker = new EvidenceTracker();
		List<KmerSupportNode> list = Stream.of(e, e2)
				.flatMap(ev -> IntStream.range(0, ev.length()).mapToObj(i -> ev.node(i)))
				.collect(Collectors.toList());
		list.stream().forEach(ksn -> tracker.track(ksn));
		Set<KmerEvidence> result = tracker.untrack(ImmutableList.of(new KmerPathSubnode(KPN(k, "AACC", 1, 1, true))));
		assertEquals(2, result.size());
	}
	@Test
	public void traverse_should_remove_both_anchor_and_nonreference_KmerEvidence_for_read_pair() {
		int k = 4;
		SAMRecord[] rp = OEA(3, 1, "4M", true);
		rp[0].setReadBases(B("AACC"));
		rp[1].setReadBases(B("ACGT"));
		NonReferenceReadPair nrrp = NonReferenceReadPair.create(rp[0], rp[1], SES(50));
		KmerEvidence e = KmerEvidence.create(k, nrrp);
		KmerEvidence e2 = KmerEvidence.createAnchor(k, nrrp, 0, nrrp.getEvidenceSource().getContext().getReference());
		EvidenceTracker tracker = new EvidenceTracker();
		List<KmerSupportNode> list = Stream.of(e, e2)
				.flatMap(ev -> IntStream.range(0, ev.length()).mapToObj(i -> ev.node(i)))
				.collect(Collectors.toList());
		list.stream().forEach(ksn -> tracker.track(ksn));
		Set<KmerEvidence> result = tracker.traverse(ImmutableList.of(new KmerPathSubnode(KPN(k, "AACC", 1, 1, true))), true);
		assertEquals(2, result.size());
	}
}
