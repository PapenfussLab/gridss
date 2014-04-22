package au.edu.wehi.socrates;

import static org.junit.Assert.*;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

import org.junit.Test;


public class LinearGenomicCoordinateTest {
	@Test
	public void testGetLinearCoordinate() {
		SAMSequenceDictionary dict = new SAMSequenceDictionary();
		dict.addSequence(new SAMSequenceRecord("contig1", 10));
		dict.addSequence(new SAMSequenceRecord("contig2", 20));
		LinearGenomicCoordinate c = new LinearGenomicCoordinate(dict);
		assertEquals(1, c.getLinearCoordinate(0,  1));
		assertEquals(11, c.getLinearCoordinate(1,  1));
	}
	@Test
	public void testGetLinearCoordinateByName() {
		SAMSequenceDictionary dict = new SAMSequenceDictionary();
		dict.addSequence(new SAMSequenceRecord("contig1", 10));
		dict.addSequence(new SAMSequenceRecord("contig2", 20));
		LinearGenomicCoordinate c = new LinearGenomicCoordinate(dict);
		assertEquals(1, c.getLinearCoordinate("contig1",  1));
		assertEquals(11, c.getLinearCoordinate("contig2",  1));
	}
	@SuppressWarnings("deprecation")
	@Test(expected=IllegalArgumentException.class)
	public void shouldRequireContigLengths() throws Exception {
		SAMSequenceDictionary dict = new SAMSequenceDictionary();
		dict.addSequence(new SAMSequenceRecord("contig1"));
		dict.addSequence(new SAMSequenceRecord("contig2"));
		LinearGenomicCoordinate c = new LinearGenomicCoordinate(dict);
		assertEquals(1, c.getLinearCoordinate(0,  1));
	}
	@Test
	public void testPadding() {
		SAMSequenceDictionary dict = new SAMSequenceDictionary();
		dict.addSequence(new SAMSequenceRecord("contig1", 10));
		dict.addSequence(new SAMSequenceRecord("contig2", 20));
		LinearGenomicCoordinate c = new LinearGenomicCoordinate(dict, 2);
		// start padding
		assertEquals(3, c.getLinearCoordinate("contig1",  1));
		// padding from start + between chromosomes = 11 + 2 + 2
		assertEquals(15, c.getLinearCoordinate("contig2",  1));
	}
}
