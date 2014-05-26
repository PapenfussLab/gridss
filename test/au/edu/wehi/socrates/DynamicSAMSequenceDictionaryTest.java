package au.edu.wehi.socrates;

import static org.junit.Assert.*;
import htsjdk.samtools.SAMSequenceRecord;

import org.junit.Test;

public class DynamicSAMSequenceDictionaryTest extends TestHelper {
	@Test
	public void getSequenceIndex_should_add() {
		DynamicSAMSequenceDictionary d = new DynamicSAMSequenceDictionary(getContext().getDictionary());
		int size = d.size();
		int index = d.getSequenceIndex("newContig");
		assertEquals(size, index);
		assertEquals(size + 1, d.size());
	}
	@Test
	public void getSequence_should_add() {
		DynamicSAMSequenceDictionary d = new DynamicSAMSequenceDictionary(getContext().getDictionary());
		SAMSequenceRecord seq = d.getSequence("newContig");
		assertEquals("newContig", seq.getSequenceName());
	}
}
