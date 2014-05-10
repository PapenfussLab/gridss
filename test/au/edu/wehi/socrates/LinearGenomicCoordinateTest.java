package au.edu.wehi.socrates;

import static org.junit.Assert.assertEquals;
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
	@Test
	public void GetStartLinearCoordinate() {
		SAMSequenceDictionary dict = new SAMSequenceDictionary();
		dict.addSequence(new SAMSequenceRecord("contig1", 10));
		dict.addSequence(new SAMSequenceRecord("contig2", 20));
		LinearGenomicCoordinate c = new LinearGenomicCoordinate(dict);
		assertEquals(2, c.getStartLinearCoordinate(new BreakpointLocation(0, BreakpointDirection.Forward, 2, 5, 0)));
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
	@Test
	public void getReferenceIndex_should_round_trip() {
		SAMSequenceDictionary dict = new SAMSequenceDictionary();
		dict.addSequence(new SAMSequenceRecord("contig1", 10));
		dict.addSequence(new SAMSequenceRecord("contig2", 20));
		dict.addSequence(new SAMSequenceRecord("contig3", 30));
		LinearGenomicCoordinate c = new LinearGenomicCoordinate(dict, 1);
		
		assertEquals(0, c.getReferenceIndex(c.getLinearCoordinate(0, 1)));
		assertEquals(1, c.getReferenceIndex(c.getLinearCoordinate(1, 1)));
	}
	@Test
	public void getReferencePosition_should_round_trip() {
		SAMSequenceDictionary dict = new SAMSequenceDictionary();
		dict.addSequence(new SAMSequenceRecord("contig1", 10));
		dict.addSequence(new SAMSequenceRecord("contig2", 20));
		dict.addSequence(new SAMSequenceRecord("contig3", 30));
		LinearGenomicCoordinate c = new LinearGenomicCoordinate(dict, 1);
		
		assertEquals(1, c.getReferencePosition(c.getLinearCoordinate(0, 1)));
		assertEquals(2, c.getReferencePosition(c.getLinearCoordinate(1, 2)));
		assertEquals(3, c.getReferencePosition(c.getLinearCoordinate(0, 3)));
		assertEquals(4, c.getReferencePosition(c.getLinearCoordinate(1, 4)));
		assertEquals(5, c.getReferencePosition(c.getLinearCoordinate(2, 5)));
		assertEquals(6, c.getReferencePosition(c.getLinearCoordinate(1, 6)));
	}
}
