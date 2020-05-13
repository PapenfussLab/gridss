package au.edu.wehi.idsv.alignment;

import au.edu.wehi.idsv.*;
import au.edu.wehi.idsv.picard.InMemoryReferenceSequenceFile;
import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.picard.SynchronousReferenceLookupAdapter;
import au.edu.wehi.idsv.picard.TwoBitBufferedReferenceSequenceFile;
import au.edu.wehi.idsv.util.SequenceUtil;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import org.junit.Assert;
import org.junit.Test;
import org.junit.experimental.categories.Category;

import java.io.File;
import java.io.FileNotFoundException;

import static org.junit.Assert.assertEquals;


public class BreakpointHomologyTest extends TestHelper {
	@Test(expected=IllegalArgumentException.class)
	public void should_require_exact_local_breakpoint() {
		InMemoryReferenceSequenceFile ref = new InMemoryReferenceSequenceFile(
				new String[] { "0", "1", },
				new byte[][] { B("CCCAATGGGCCC"),
							   B("TTTAATGGGAAA"), });
								//     >
								//      <
		BreakpointHomology.calculate(ref, new BreakpointSummary(0, FWD, 6, 6, 7, 1, BWD, 7, 7, 7), "", 100, 0);
	}
	@Test(expected=IllegalArgumentException.class)
	public void should_require_exact_remote_breakpoint() {
		InMemoryReferenceSequenceFile ref = new InMemoryReferenceSequenceFile(
				new String[] { "0", "1", },
				new byte[][] { B("CCCAATGGGCCC"),
							   B("TTTAATGGGAAA"), });
								//     >
								//      <
		BreakpointHomology.calculate(ref, new BreakpointSummary(0, FWD, 6, 6, 6, 1, BWD, 7, 7, 8), "", 100, 0);
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
		BreakpointHomology bh = BreakpointHomology.calculate(ref, new BreakpointSummary(0, FWD, 6, 1, BWD, 7), "", 100, 0);
		assertEquals(6, bh.getRemoteHomologyLength() + bh.getLocalHomologyLength());
		assertEquals(3, bh.getLocalHomologyLength());
		assertEquals(3, bh.getRemoteHomologyLength());
	}
	@Test
	public void no_homology() {
		InMemoryReferenceSequenceFile ref = new InMemoryReferenceSequenceFile(
				new String[] { "0", "1", },
				new byte[][] { B("CCCCCCCCCCCC"),
							   B("TTTTTTTTTTTTTTTTTTTTTTTT"), });
								//     >
								//      <
		BreakpointHomology bh = BreakpointHomology.calculate(ref, new BreakpointSummary(0, FWD, 6, 1, BWD, 7), "", 100, 10);
		assertEquals(0, bh.getLocalHomologyLength());
		assertEquals(0, bh.getRemoteHomologyLength());
	}
	@Test
	public void should_report_homology_up_to_max_distance_away() {
		InMemoryReferenceSequenceFile ref = new InMemoryReferenceSequenceFile(
				new String[] { "0", "1", },
				new byte[][] { B("CCCAATGGGCCCTTTTTTTTTTTTTTTTTTT"),
							   B("TTTAATGGGAAAGGGGGGGGGGGGGGGGGGG"), });
								//     >
								//      <
		BreakpointHomology bh = BreakpointHomology.calculate(ref, new BreakpointSummary(0, FWD, 6, 1, BWD, 7), "", 1, 5);
		assertEquals(1, bh.getLocalHomologyLength());
		assertEquals(1, bh.getRemoteHomologyLength());
	}
	@Test
	public void should_reduce_window_size_for_small_events() {
		BreakpointHomology bh = BreakpointHomology.calculate(SMALL_FA, new BreakpointSummary(0, FWD, 100, 0, BWD, 103), "", 100, 0);
		assertEquals(2, bh.getLocalHomologyLength());
		assertEquals(2, bh.getRemoteHomologyLength());
	}
	@Test
	public void should_incorporate_inserted_sequence() {
		InMemoryReferenceSequenceFile ref = new InMemoryReferenceSequenceFile(
				new String[] { "0", "1", },
				new byte[][] { B("CCCAAAATTTGGGAAAAAATTTTTTTTTTTTTTTTTTT"),
							   B("TTTAAAATTTGGGAAAAAAGGGGGGGGGGGGGGGGGGG"), });
		// CCCAAAATTT GGGAAAAAATTTTTTTTTTT TTTTTTTT"),
		// TTTAAAATTT GGGAAAAAAGGGGGGGGGGG GGGGGGGG"
		// MMMMMMMMMMIMMMMMMMMSSSSSSSSSSSS fwd
		// SSSMMMMMMMIMMMMMMMMMMMMMMMMMMMM bwd
		// 1234567890 1234567890
		//          >G
		//            <
		// 
		BreakpointHomology bh = BreakpointHomology.calculate(ref, new BreakpointSummary(0, FWD, 10, 1, BWD, 11), "c", 20, 10);
		assertEquals(7, bh.getLocalHomologyLength());
		assertEquals(9, bh.getRemoteHomologyLength());
	}
	@Test
	public void should_not_exceed_contig_bounds() {
		InMemoryReferenceSequenceFile underlying = new InMemoryReferenceSequenceFile(
				new String[] { "0", "1", },
				new byte[][] { B("CCCCCCCCCCC"),
							   B("AAAAAAAAAAA"), });
		//BreakpointHomology bh0 = BreakpointHomology.calculate(ref, new BreakpointSummary(0, FWD, 1, 1, BWD, 20), "", 50, 50);
		TwoBitBufferedReferenceSequenceFile ref = new TwoBitBufferedReferenceSequenceFile(underlying);
		for (int b1pos = 1; b1pos <= ref.getSequenceDictionary().getSequences().get(0).getSequenceLength(); b1pos++) {
			for (int b2pos = 1; b2pos <= ref.getSequenceDictionary().getSequences().get(0).getSequenceLength(); b2pos++) {
				for (BreakendDirection b1dir : BreakendDirection.values()) {
					for (BreakendDirection b2dir : BreakendDirection.values()) {
						for (int maxBreakendLength = 1; maxBreakendLength < ref.getSequenceDictionary().getSequences().get(0).getSequenceLength() + 2; maxBreakendLength++) {
							for (int margin = 0; margin < ref.getSequenceDictionary().getSequences().get(0).getSequenceLength() + 2; margin++) {
								BreakpointHomology bh = BreakpointHomology.calculate(ref, new BreakpointSummary(0, b1dir, b1pos, 1, b2dir, b2pos), "", maxBreakendLength, margin);
								assertEquals(0, bh.getLocalHomologyLength());
								assertEquals(0, bh.getRemoteHomologyLength());
							}
						}
					}
				}
			}
		}
	}
	@Test
	public void should_not_report_homology_for_small_dup() {
		InMemoryReferenceSequenceFile ref = new InMemoryReferenceSequenceFile(
				new String[] { "0", "1", },
				new byte[][] {
						//    0        1         2         3         4         5         6
						//    123456789012345678901234567890123456789012345678901234567890
						B("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACC"),
						B("CACCGCCTATTCGAACGGGCGAATCTACCTAGGTCGCTCAGAACCGGCACCCTTAA"), });

		BreakpointHomology bh = BreakpointHomology.calculate(ref, new BreakpointSummary(0, BWD, 10, 0, FWD, 20), "", 20, 10);
		assertEquals(0, bh.getLocalHomologyLength());
		assertEquals(0, bh.getRemoteHomologyLength());
	}
	@Test
	public void should_handle_inserted_sequence_no_homology() {
		InMemoryReferenceSequenceFile ref = new InMemoryReferenceSequenceFile(
				new String[] { "0", "1", },
				new byte[][] {
						//    0        1         2         3         4         5         6
						//    123456789012345678901234567890123456789012345678901234567890
						B("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACC"),
						B("CACCGCCTATTCGAACGGGCGAATCTACCTAGGTCGCTCAGAACCGGCACCCTTAA"), });

		BreakpointHomology bh = BreakpointHomology.calculate(ref, new BreakpointSummary(0, FWD, 10, 1, BWD, 10), "GGGGG", 20, 10);
		assertEquals(0, bh.getLocalHomologyLength());
		assertEquals(0, bh.getRemoteHomologyLength());
	}
	@Test
	public void should_handle_inserted_sequence_with_homology() {
		InMemoryReferenceSequenceFile ref = new InMemoryReferenceSequenceFile(
				new String[] { "0", "1", },
				new byte[][] {
						//    0        1         2         3         4         5         6
						//    123456789012345678901234567890123456789012345678901234567890
						B("AGAACCGGCACCCTTAACATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACC"),
						B("AGAACCGGCACCCTTAACACCGCCTACAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACC"), });

		BreakpointSummary bp = new BreakpointSummary(0, FWD, 27, 1, BWD, 28);
		BreakpointHomology bh = BreakpointHomology.calculate(ref, bp, "GGGGG", 20, 5);
		assertEquals(0, bh.getLocalHomologyLength());
		assertEquals(20, bh.getRemoteHomologyLength());

		bh = BreakpointHomology.calculate(ref, bp.remoteBreakpoint(), "GGGGG", 20, 5);
		assertEquals(0, bh.getRemoteHomologyLength());
		assertEquals(20, bh.getLocalHomologyLength());
	}
	@Test
	@Category(Hg19Tests.class)
	public void issue344_regression_ihompos_should_be_symmetrical() throws FileNotFoundException {
		InMemoryReferenceSequenceFile zref = new InMemoryReferenceSequenceFile(
				new String[] { "0", "1", },
				new byte[][] {
						//    0        1         2         3         4         5         6
						//    123456789012345678901234567890123456789012345678901234567890
						B("AGAACCGGCACCCTTAACATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACC"),
						B("AGAACCGGCACCCTTAACACCGCCTACAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACC"), });

		IndexedFastaSequenceFile ref = new IndexedFastaSequenceFile(Hg19Tests.findBroadHg19Reference());
		ReferenceLookup hg19ref = new SynchronousReferenceLookupAdapter(ref);
		SAMSequenceDictionary dict = ref.getSequenceDictionary();
		BreakpointSummary bs = new BreakpointSummary(dict.getSequenceIndex("1"), FWD, 145396172, dict.getSequenceIndex("1"), BWD, 145396592);
		BreakpointHomology bhLeft = BreakpointHomology.calculate(hg19ref, bs, "TC", 20, 5);
		BreakpointHomology bhRight = BreakpointHomology.calculate(hg19ref, bs.remoteBreakpoint(), "TC", 20, 5);
		Assert.assertEquals(bhLeft.getLocalHomologyLength(), bhRight.getRemoteHomologyLength());
		Assert.assertEquals(bhLeft.getRemoteHomologyLength(), bhRight.getLocalHomologyLength());
	}

	/**
	 * Technically we should report the full homology length but that's tricky to calculate
	 */
	@Test
	public void should_report_single_instance_homology_for_small_event_when_ref_is_homologous() {
		InMemoryReferenceSequenceFile ref = new InMemoryReferenceSequenceFile(
				new String[] { "0",},
				new byte[][] {
						//    0        1         2         3         4         5         6
						//    123456789012345678901234567890123456789012345678901234567890
						B("CATTAATCGCAAGAGCGGGCAAGAGCGGGCAAGAGCGGGCAAGAGCGGGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACC")});

		BreakpointHomology bh = BreakpointHomology.calculate(ref, new BreakpointSummary(0, FWD, 30, 0, BWD, 41), "", 20, 10);
		assertEquals(10, bh.getLocalHomologyLength());
		assertEquals(10, bh.getRemoteHomologyLength());
	}

	@Test
	public void should_report_insertion_homology() {
		BreakpointSummary bp = new BreakpointSummary(0, FWD, 100, 0, BWD, 500);
		BreakpointHomology bh = BreakpointHomology.calculate(SMALL_FA, bp, "AAAAA", 10, 50);
		assertEquals(10, bh.getLocalHomologyLength());
		assertEquals(10, bh.getRemoteHomologyLength());
		bh = BreakpointHomology.calculate(SMALL_FA, bp.remoteBreakpoint(), "AAAAA", 10, 50);
		assertEquals(10, bh.getLocalHomologyLength());
		assertEquals(10, bh.getRemoteHomologyLength());
	}
}
