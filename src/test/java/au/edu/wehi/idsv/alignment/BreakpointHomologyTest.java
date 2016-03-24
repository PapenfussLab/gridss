package au.edu.wehi.idsv.alignment;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import au.edu.wehi.idsv.BreakpointSummary;
import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.picard.InMemoryReferenceSequenceFile;


public class BreakpointHomologyTest extends TestHelper {
	@Test(expected=IllegalArgumentException.class)
	public void should_require_exact_local_breakpoint() {
		InMemoryReferenceSequenceFile ref = new InMemoryReferenceSequenceFile(
				new String[] { "0", "1", },
				new byte[][] { B("CCCAATGGGCCC"),
							   B("TTTAATGGGAAA"), });
								//     >
								//      <
		BreakpointHomology.calculate(ref, new BreakpointSummary(0, FWD, 6, 7, 1, BWD, 7, 7), "", 100, 0);
	}
	@Test(expected=IllegalArgumentException.class)
	public void should_require_exact_remote_breakpoint() {
		InMemoryReferenceSequenceFile ref = new InMemoryReferenceSequenceFile(
				new String[] { "0", "1", },
				new byte[][] { B("CCCAATGGGCCC"),
							   B("TTTAATGGGAAA"), });
								//     >
								//      <
		BreakpointHomology.calculate(ref, new BreakpointSummary(0, FWD, 6, 6, 1, BWD, 7, 8), "", 100, 0);
	}
	@Test
	public void should_report_homology_on_both_sides() {
		InMemoryReferenceSequenceFile ref = new InMemoryReferenceSequenceFile(
				new String[] { "0", "1", },
				new byte[][] { B("CCCAATGGGCCC"),
							   B("TTTAATGGGAAA"), });
								//123456789012
								//     >
								//      <
		BreakpointHomology bh = BreakpointHomology.calculate(ref, new BreakpointSummary(0, FWD, 6, 6, 1, BWD, 7, 7), "", 100, 0);
		assertEquals(6, bh.getUpperBound() + bh.getLowerBound());
		assertEquals(3, bh.getLowerBound());
		assertEquals(3, bh.getUpperBound());
	}
	@Test
	public void no_homology() {
		InMemoryReferenceSequenceFile ref = new InMemoryReferenceSequenceFile(
				new String[] { "0", "1", },
				new byte[][] { B("CCCCCCCCCCCC"),
							   B("TTTTTTTTTTTT"), });
								//     >
								//      <
		BreakpointHomology bh = BreakpointHomology.calculate(ref, new BreakpointSummary(0, FWD, 6, 6, 1, BWD, 7, 7), "", 100, 0);
		assertEquals(0, bh.getLowerBound());
		assertEquals(0, bh.getUpperBound());
	}
	@Test
	public void should_report_homology_up_to_max_distance_away() {
		InMemoryReferenceSequenceFile ref = new InMemoryReferenceSequenceFile(
				new String[] { "0", "1", },
				new byte[][] { B("CCCAATGGGCCC"),
							   B("TTTAATGGGAAA"), });
								//     >
								//      <
		BreakpointHomology bh = BreakpointHomology.calculate(ref, new BreakpointSummary(0, FWD, 6, 6, 1, BWD, 7, 7), "", 1, 0);
		assertEquals(1, bh.getLowerBound());
		assertEquals(1, bh.getUpperBound());
	}
	@Test
	public void should_reduce_window_size_for_small_events() {
		BreakpointHomology bh = BreakpointHomology.calculate(SMALL_FA, new BreakpointSummary(0, FWD, 100, 100, 0, BWD, 103, 103), "", 100, 0);
		assertEquals(2, bh.getLowerBound());
		assertEquals(2, bh.getUpperBound());
	}
	@Test
	public void should_incorporate_inserted_sequence() {
		InMemoryReferenceSequenceFile ref = new InMemoryReferenceSequenceFile(
				new String[] { "0", "1", },
				new byte[][] { B("CCCAATGGGCCC"),
							   B("TTTAATGGGAAA"), });
								//     >G
								//       <
		BreakpointHomology bh = BreakpointHomology.calculate(ref, new BreakpointSummary(0, FWD, 6, 6, 1, BWD, 8, 8), "G", 10, 0);
		assertEquals(3, bh.getLowerBound());
		assertEquals(3, bh.getUpperBound());
	}
}
