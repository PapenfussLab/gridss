package au.edu.wehi.idsv.sam;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;

import java.util.ArrayList;

import org.junit.Test;

import com.google.common.collect.Lists;
import com.google.common.collect.ImmutableList;


public class SAMRecordCigarCleaningIteratorTest {
	@Test
	public void should_clean_cigar() {
		SAMRecord r = new SAMRecord(null);
		r.setCigarString("10I90M10M");
		ArrayList<SAMRecord> list = Lists.newArrayList(new SAMRecordCigarCleaningIterator(ImmutableList.of(r).iterator()));
		assertEquals("10S100M", list.get(0).getCigarString());
	}
	@Test
	public void should_not_touch_correct_cigar() {
		SAMRecord r1 = new SAMRecord(null);
		SAMRecord r2 = new SAMRecord(null);
		r2.setCigarString("100M");
		Cigar c1 = r1.getCigar();
		Cigar c2 = r2.getCigar();
		ArrayList<SAMRecord> list = Lists.newArrayList(new SAMRecordCigarCleaningIterator(ImmutableList.of(r1, r2).iterator()));
		assertTrue(list.get(0).getCigar() == c1);
		assertTrue(list.get(1).getCigar() == c2);
	}
}
