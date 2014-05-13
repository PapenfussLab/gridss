package au.edu.wehi.socrates;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.util.List;

import htsjdk.samtools.SAMRecord;

import org.junit.Test;

import com.google.common.collect.Iterators;

public class SequentialNonReferenceReadPairFactoryTest extends TestHelper {
	public SequentialNonReferenceReadPairFactory getFactory(List<SAMRecord> data) {
		data = mateSorted(data);
		return new SequentialNonReferenceReadPairFactory(Iterators.peekingIterator(data.iterator()));
	}
	@Test
	public void should_match_by_read_name_and_pair_index() {
		SequentialNonReferenceReadPairFactory factory = getFactory(toList(OEA(0, 1, "100M", true)[1]));
		NonReferenceReadPair rp = factory.createNonReferenceReadPair(OEA(0, 1, "100M", true)[0], 100);
		assertNotNull(rp);
	}
	@Test
	public void should_not_match_with_read_name_mismatch() {
		SequentialNonReferenceReadPairFactory factory = getFactory(toList(OEA(0, 1, "100M", true)[1]));
		NonReferenceReadPair rp = factory.createNonReferenceReadPair(withReadName("different", OEA(0, 1, "100M", true))[0], 100);
		assertNull(rp);
	}
	@Test
	public void should_not_match_with_read_pair_index_mismatch() {
		SequentialNonReferenceReadPairFactory factory = getFactory(toList(OEA(0, 1, "100M", true)[1]));
		SAMRecord r = OEA(0, 1, "100M", true)[0];
		r.setFirstOfPairFlag(!r.getFirstOfPairFlag());
		NonReferenceReadPair rp = factory.createNonReferenceReadPair(r, 100);
		assertNull(rp);
	}
	@Test
	public void should_match_sequential_records() {
		SAMRecord[] mates = new SAMRecord[] {
			withReadName("r1", OEA(0, 1, "100M", true))[1],
			withReadName("r3", OEA(0, 1, "100M", true))[1],
			withReadName("r2", OEA(0, 1, "100M", true))[1], // equivalence class same pos
			withReadName("r4", OEA(0, 10, "100M", true))[1], // equivalence class advance pos
			withReadName("r5", OEA(1, 1, "100M", true))[1], // equivalence class down pos
			withReadName("r6", OEA(2, 1, "100M", true))[1], // equivalence class equal pos
			withReadName("r7", OEA(3, 3, "100M", true))[1], // equivalence class up pos
		};
		SequentialNonReferenceReadPairFactory factory = getFactory(toList(mates));
		assertEquals("r1", factory.createNonReferenceReadPair(withReadName("r1", OEA(0, 1, "100M", true))[0], 100).getEvidenceID());
		// r2 and r3 should be retrievable since they are 
		assertEquals("r2", factory.createNonReferenceReadPair(withReadName("r2", OEA(0, 1, "100M", true))[0], 100).getEvidenceID());
		assertEquals("r3", factory.createNonReferenceReadPair(withReadName("r3", OEA(0, 1, "100M", true))[0], 100).getEvidenceID());
		assertEquals("r4", factory.createNonReferenceReadPair(withReadName("r4", OEA(0, 10, "100M", true))[0], 100).getEvidenceID());
		assertEquals("r5", factory.createNonReferenceReadPair(withReadName("r5", OEA(1, 1, "100M", true))[0], 100).getEvidenceID());
		assertEquals("r6", factory.createNonReferenceReadPair(withReadName("r6", OEA(2, 1, "100M", true))[0], 100).getEvidenceID());
		assertEquals("r7", factory.createNonReferenceReadPair(withReadName("r7", OEA(3, 3, "100M", true))[0], 100).getEvidenceID());
	}
	@Test(expected=IllegalStateException.class)
	public void should_fail_during_non_sequential_traversal() {
		SAMRecord[] mates = new SAMRecord[] {
			withReadName("r1", OEA(0, 2, "100M", true))[1],
			withReadName("r2", OEA(0, 1, "100M", true))[1],
		};
		SequentialNonReferenceReadPairFactory factory = getFactory(toList(mates));
		assertEquals("r1", factory.createNonReferenceReadPair(withReadName("r1", OEA(0, 2, "100M", true))[0], 100).getEvidenceID()); 
		assertNull(factory.createNonReferenceReadPair(withReadName("r2", OEA(0, 1, "100M", true))[0], 100));
	}
}
