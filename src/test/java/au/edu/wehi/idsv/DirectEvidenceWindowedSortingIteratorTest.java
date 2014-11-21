package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

import java.util.ArrayList;
import java.util.Collections;

import org.junit.Test;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;


public class DirectEvidenceWindowedSortingIteratorTest extends TestHelper {
	@Test
	public void should_sort_by_natural_ordering() {
		ImmutableList<DirectedEvidence> list = ImmutableList.of(
				E(0, 1, FWD),
				E(0, 6, FWD),
				E(0, 5, FWD),
				E(0, 4, FWD),
				E(0, 7, FWD),
				E(0, 7, 7, BWD),
				E(0, 7, 8, BWD),
				E(0, 7, 6, BWD)
				);
		DirectEvidenceWindowedSortingIterator<DirectedEvidence> it = new DirectEvidenceWindowedSortingIterator<DirectedEvidence>(getContext(), 5, list.iterator());
		ArrayList<DirectedEvidence> expected = Lists.newArrayList(list);
		Collections.sort(expected, DirectedEvidenceOrder.ByNatural);
		ArrayList<DirectedEvidence> result = Lists.newArrayList(it);
		assertEquals(expected, result);
	}
	@Test
	public void should_not_be_able_to_sort_unordered_outside_of_window() {
		ImmutableList<DirectedEvidence> list = ImmutableList.of(
				E(0, 10, FWD),
				E(0, 20, FWD),
				E(0, 30, FWD),
				E(0, 50, FWD),
				E(0, 60, FWD),
				E(0, 40, FWD)
				);
		DirectEvidenceWindowedSortingIterator<DirectedEvidence> it = new DirectEvidenceWindowedSortingIterator<DirectedEvidence>(getContext(), 5, list.iterator());
		ArrayList<DirectedEvidence> expected = Lists.newArrayList(list);
		Collections.sort(expected, DirectedEvidenceOrder.ByNatural);
		ArrayList<DirectedEvidence> result = Lists.newArrayList(it);
		assertFalse(expected.equals(result));
	}
}
