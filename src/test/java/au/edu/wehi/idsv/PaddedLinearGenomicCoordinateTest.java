package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;


public class PaddedLinearGenomicCoordinateTest {
	@Test
	public void testGetLinearCoordinate() {
		SAMSequenceDictionary dict = new SAMSequenceDictionary();
		dict.addSequence(new SAMSequenceRecord("contig1", 10));
		dict.addSequence(new SAMSequenceRecord("contig2", 20));
		LinearGenomicCoordinate c = new PaddedLinearGenomicCoordinate(dict);
		assertEquals(1, c.getLinearCoordinate(0,  1));
		assertEquals(11, c.getLinearCoordinate(1,  1));
	}
	@Test
	public void testGetLinearCoordinateByName() {
		SAMSequenceDictionary dict = new SAMSequenceDictionary();
		dict.addSequence(new SAMSequenceRecord("contig1", 10));
		dict.addSequence(new SAMSequenceRecord("contig2", 20));
		LinearGenomicCoordinate c = new PaddedLinearGenomicCoordinate(dict);
		assertEquals(1, c.getLinearCoordinate("contig1",  1));
		assertEquals(11, c.getLinearCoordinate("contig2",  1));
	}
	@Test
	public void fixed_padding_should_offset_by_fixed_amount() {
		SAMSequenceDictionary dict = new SAMSequenceDictionary();
		dict.addSequence(new SAMSequenceRecord("contig1", 10));
		dict.addSequence(new SAMSequenceRecord("contig2", 20));
		LinearGenomicCoordinate c = new PaddedLinearGenomicCoordinate(dict, 100000, true);
		assertEquals(100001, c.getLinearCoordinate("contig1",  1));
		assertEquals(200001, c.getLinearCoordinate("contig2",  1));
	}
	@Test(expected=IllegalArgumentException.class)
	public void fixed_padding_should_not_be_smaller_than_max_contig_length() {
		SAMSequenceDictionary dict = new SAMSequenceDictionary();
		dict.addSequence(new SAMSequenceRecord("contig1", 10));
		dict.addSequence(new SAMSequenceRecord("contig2", 20));
		new PaddedLinearGenomicCoordinate(dict, 19, true);
	}
	@Test
	public void GetStartLinearCoordinate() {
		SAMSequenceDictionary dict = new SAMSequenceDictionary();
		dict.addSequence(new SAMSequenceRecord("contig1", 10));
		dict.addSequence(new SAMSequenceRecord("contig2", 20));
		LinearGenomicCoordinate c = new PaddedLinearGenomicCoordinate(dict);
		assertEquals(2, c.getStartLinearCoordinate(new BreakendSummary(0, BreakendDirection.Forward, 2, 5)));
	}
	@Test
	public void GetEndLinearCoordinate() {
		SAMSequenceDictionary dict = new SAMSequenceDictionary();
		dict.addSequence(new SAMSequenceRecord("contig1", 10));
		dict.addSequence(new SAMSequenceRecord("contig2", 20));
		LinearGenomicCoordinate c = new PaddedLinearGenomicCoordinate(dict);
		assertEquals(5, c.getEndLinearCoordinate(new BreakendSummary(0, BreakendDirection.Forward, 2, 5)));
	}
	@SuppressWarnings("deprecation")
	@Test(expected=IllegalArgumentException.class)
	public void shouldRequireContigLengths() throws Exception {
		SAMSequenceDictionary dict = new SAMSequenceDictionary();
		dict.addSequence(new SAMSequenceRecord("contig1"));
		dict.addSequence(new SAMSequenceRecord("contig2"));
		LinearGenomicCoordinate c = new PaddedLinearGenomicCoordinate(dict);
		assertEquals(1, c.getLinearCoordinate(0,  1));
	}
	@Test
	public void testPadding() {
		SAMSequenceDictionary dict = new SAMSequenceDictionary();
		dict.addSequence(new SAMSequenceRecord("contig1", 10));
		dict.addSequence(new SAMSequenceRecord("contig2", 20));
		LinearGenomicCoordinate c = new PaddedLinearGenomicCoordinate(dict, 2);
		// start padding
		assertEquals(3, c.getLinearCoordinate("contig1",  1));
		// padding from start + between chromosomes = 2 + 11 + 2 + 2
		assertEquals(17, c.getLinearCoordinate("contig2",  1));
	}
	@Test
	public void getReferenceIndex_should_round_trip() {
		SAMSequenceDictionary dict = new SAMSequenceDictionary();
		dict.addSequence(new SAMSequenceRecord("contig1", 10));
		dict.addSequence(new SAMSequenceRecord("contig2", 20));
		dict.addSequence(new SAMSequenceRecord("contig3", 30));
		LinearGenomicCoordinate c = new PaddedLinearGenomicCoordinate(dict, 1);
		
		assertEquals(0, c.getReferenceIndex(c.getLinearCoordinate(0, 1)));
		assertEquals(1, c.getReferenceIndex(c.getLinearCoordinate(1, 1)));
	}
	@Test
	public void getReferencePosition_should_round_trip() {
		SAMSequenceDictionary dict = new SAMSequenceDictionary();
		dict.addSequence(new SAMSequenceRecord("contig1", 10));
		dict.addSequence(new SAMSequenceRecord("contig2", 20));
		dict.addSequence(new SAMSequenceRecord("contig3", 30));
		for (int bufferSize = 0; bufferSize < 3; bufferSize++) {
			LinearGenomicCoordinate c = new PaddedLinearGenomicCoordinate(dict, bufferSize);
			assertEquals(1, c.getReferencePosition(c.getLinearCoordinate(0, 1)));
			assertEquals(2, c.getReferencePosition(c.getLinearCoordinate(1, 2)));
			assertEquals(3, c.getReferencePosition(c.getLinearCoordinate(0, 3)));
			assertEquals(4, c.getReferencePosition(c.getLinearCoordinate(1, 4)));
			assertEquals(5, c.getReferencePosition(c.getLinearCoordinate(2, 5)));
			assertEquals(6, c.getReferencePosition(c.getLinearCoordinate(1, 6)));
		}
	}
	@Test
	public void getReferencePosition_should_round_trip_all_chr_positions() {
		SAMSequenceDictionary dict = new SAMSequenceDictionary();
		dict.addSequence(new SAMSequenceRecord("contig1", 10));
		dict.addSequence(new SAMSequenceRecord("contig2", 10));
		dict.addSequence(new SAMSequenceRecord("contig3", 10));
		for (int bufferSize = 0; bufferSize < 3; bufferSize++) {
			LinearGenomicCoordinate c = new PaddedLinearGenomicCoordinate(dict, bufferSize);
			for (int i = 1; i <= 10; i++) {
				assertEquals(i, c.getReferencePosition(c.getLinearCoordinate(0, i)));
				assertEquals(i, c.getReferencePosition(c.getLinearCoordinate(1, i)));
				assertEquals(i, c.getReferencePosition(c.getLinearCoordinate(2, i)));
			}
		}
	}
	@Test
	public void should_round_trip_before_and_after_chr_positions() {
		SAMSequenceDictionary dict = new SAMSequenceDictionary();
		dict.addSequence(new SAMSequenceRecord("contig1", 10));
		dict.addSequence(new SAMSequenceRecord("contig2", 10));
		dict.addSequence(new SAMSequenceRecord("contig3", 10));
		for (int bufferSize = 0; bufferSize < 16; bufferSize++) {
			LinearGenomicCoordinate c = new PaddedLinearGenomicCoordinate(dict, bufferSize);
			for (int i = 1 - bufferSize; i <= 10 + bufferSize; i++) {
				for (int j = 0; j < 3; j++) {
					assertEquals(j, c.getReferenceIndex(c.getLinearCoordinate(j, i)));
					assertEquals(i, c.getReferencePosition(c.getLinearCoordinate(j, i)));
				}
			}
		}
	}
	@Test
	public void should_map_to_closest_chr_fixed_width() {
		SAMSequenceDictionary dict = new SAMSequenceDictionary();
		dict.addSequence(new SAMSequenceRecord("contig0", 1));
		dict.addSequence(new SAMSequenceRecord("contig1", 2));
		dict.addSequence(new SAMSequenceRecord("contig2", 3));
		//            1         2         3
		//  0123456789012345678901234567890
		//  ++++0+++11++222+++++
		// --00000*1111222222222-
		//  aaaabbbbccccddddeeee   width blocks
		LinearGenomicCoordinate c = new PaddedLinearGenomicCoordinate(dict, 4, true);
		for (int i = -10; i < 30; i++) {
			int refIndex = c.getReferenceIndex(i);
			if (i <= 0) assertEquals(-1, refIndex);
			else if (i < 7) assertEquals(0, refIndex);
			else if (i < 11) assertEquals(1, refIndex);
			else if (i < 20) assertEquals(2, refIndex);
			else assertEquals(-1, refIndex);
		}
	}
	@Test
	public void should_map_to_closest_chr_padded() {
		SAMSequenceDictionary dict = new SAMSequenceDictionary();
		dict.addSequence(new SAMSequenceRecord("contig0", 1));
		dict.addSequence(new SAMSequenceRecord("contig1", 2));
		dict.addSequence(new SAMSequenceRecord("contig2", 3));
		//            1         2         3
		//  01234567890123456789012345678901234567890
		//   ++++0++++++++11++++++++222++++
		// --000000000111111111122222222222-
		//  aaaa aaaabbbb  bbbbcccc   cccc
		LinearGenomicCoordinate c = new PaddedLinearGenomicCoordinate(dict, 4, false);
		for (int i = -10; i < 40; i++) {
			int refIndex = c.getReferenceIndex(i);
			if (i <= 0) assertEquals(-1, refIndex);
			else if (i <= 9) assertEquals(0, refIndex);
			else if (i <= 19) assertEquals(1, refIndex);
			else if (i <= 30) assertEquals(2, refIndex);
			else assertEquals(-1, refIndex);
		}
	}
	@Test
	public void getReferencePosition_should_be_negative_one_outside_of_padding() {
		SAMSequenceDictionary dict = new SAMSequenceDictionary();
		dict.addSequence(new SAMSequenceRecord("contig1", 10));
		dict.addSequence(new SAMSequenceRecord("contig2", 10));
		LinearGenomicCoordinate c = new PaddedLinearGenomicCoordinate(dict, 1);
		assertEquals(-1, c.getReferenceIndex(0));
		assertEquals(-1, c.getReferenceIndex(1+10+1 +1+10+1 + 2));
	}
}
