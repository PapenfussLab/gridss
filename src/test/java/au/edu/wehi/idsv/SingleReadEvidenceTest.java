package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Collection;
import java.util.List;
import java.util.Map.Entry;
import java.util.stream.Collectors;

import com.google.common.collect.Range;
import org.junit.Ignore;
import org.junit.Test;
import org.junit.experimental.categories.Category;

import com.google.common.collect.ImmutableListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimaps;

import au.edu.wehi.idsv.picard.InMemoryReferenceSequenceFile;
import au.edu.wehi.idsv.picard.SynchronousReferenceLookupAdapter;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.SequenceUtil;

public class SingleReadEvidenceTest extends TestHelper {
	@Test
	public void should_not_create_indel_for_XNX_placeholder() {
		for (SAMRecord r : new SAMRecord[] {
				Read(0, 1, "1X10S"),
				Read(0, 1, "2X10S"),
				Read(0, 1, "1X1N1X10S"),
				Read(0, 1, "1X100N1X10S"),
		}) {
			List<SingleReadEvidence> list = SingleReadEvidence.createEvidence(SES(), 0, r);
			assertEquals(1, list.size());
			assertTrue(list.get(0) instanceof SoftClipEvidence);
		}
	}
	@Test
	public void should_not_calculate_homology_if_untemplated_sequence_exists() {
		SAMRecord realign = Read(0, 40, "10S10M");
		realign.setAttribute("SA", "polyA,10,+,9M1S10S,0,0");
		SplitReadEvidence e = SplitReadEvidence.create(SES(), realign).get(0);
		assertEquals("", e.getHomologySequence());
		assertEquals(0, e.getHomologyAnchoredBaseCount());
	}
	@Test
	public void should_not_calculate_homology_if_inexact() {
		SAMRecord realign = Read(0, 1, "1X10S");
		realign.setAttribute("SA", "polyA,10,+,1S10M,0,0");
		SplitReadEvidence e = SplitReadEvidence.create(SES(), realign).get(0);
		assertFalse(e.isBreakendExact());
		assertEquals("", e.getHomologySequence());
		assertEquals(0, e.getHomologyAnchoredBaseCount());
	}
	public void microhomology_match_fb() {
		// 0        1         2         3         4         5         6         7         
		// 1234567890123456789012345678901234567890123456789012345678901234567890123456789
		// ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA
		//     >>>>    <<<<     
		SAMRecord r = withAttr("SA", "polyACGT,13,+,4S4M,0,0", withSequence("ACGTACGN", Read(1, 5, "4M4S")))[0];
		SplitReadEvidence e = (SplitReadEvidence)asAssemblyEvidence(r);
		assertEquals("ACGTACG", e.getHomologySequence());
		assertEquals(4, e.getHomologyAnchoredBaseCount());
		assertEquals(8, e.getBreakendSummary().nominal);
		assertEquals(13, e.getBreakendSummary().nominal2);
		assertEquals(4, e.getBreakendSummary().start);
		assertEquals(10, e.getBreakendSummary().end);
		assertEquals(16, e.getBreakendSummary().start2);
		assertEquals(9, e.getBreakendSummary().end2);
	}
	public void microhomology_match_consider_strand() {
		// 0        1         2         3         4         5         6         7         
		// 1234567890123456789012345678901234567890123456789012345678901234567890123456789
		// ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA
		//     >>>>    <<<<     
		SAMRecord r = withAttr("SA", "polyACGT,13,-,4S4M,0,0", onNegative(withSequence("ACGTACGN", Read(1, 5, "4M4S"))))[0];
		SplitReadEvidence e = (SplitReadEvidence)asAssemblyEvidence(r);
		assertEquals("ACGTACG", e.getHomologySequence());
		assertEquals(4, e.getHomologyAnchoredBaseCount());
		assertEquals(8, e.getBreakendSummary().nominal);
		assertEquals(13, e.getBreakendSummary().nominal2);
		assertEquals(4, e.getBreakendSummary().start);
		assertEquals(10, e.getBreakendSummary().end);
		assertEquals(16, e.getBreakendSummary().start2);
		assertEquals(9, e.getBreakendSummary().end2);
	}
	public void microhomology_match_bf() {
		// 0        1         2         3         4         5         6         7         
		// 1234567890123456789012345678901234567890123456789012345678901234567890123456789
		// ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA
		//     >>>>    <<<<
		SAMRecord r = withAttr("SA", "polyACGT,5,+,4M4S,0,0", withSequence("ACGTACGN", Read(1, 13, "4S4M")))[0];
		SplitReadEvidence e = (SplitReadEvidence)asAssemblyEvidence(r);
		assertEquals("ACGTACG", e.getHomologySequence());
		assertEquals(3, e.getHomologyAnchoredBaseCount());
		assertEquals(8, e.getBreakendSummary().nominal2);
		assertEquals(13, e.getBreakendSummary().nominal);
		assertEquals(4, e.getBreakendSummary().start2);
		assertEquals(10, e.getBreakendSummary().end2);
		assertEquals(16, e.getBreakendSummary().start);
		assertEquals(9, e.getBreakendSummary().end);
	}
	public void microhomology_match_ff() {
		//      1         2         3         4
		// 1234567890123456789012345678901234567890123456789
		// CATTAATCGCAATAAAAATGTTCTTTTTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAA
		//            >>>>      >>>>>
		String contig = "CATTAATCGCAATAAAAATGTTCTTTTTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAA";
		InMemoryReferenceSequenceFile ref = new InMemoryReferenceSequenceFile(new String[] { "Contig" }, new byte[][] { B(contig) });
		SAMRecord r = new SAMRecord(new SAMFileHeader());
		r.getHeader().setSequenceDictionary(getSequenceDictionary());
		r.setReferenceIndex(0);
		r.setCigarString("4M5S");
		r.setAlignmentStart(1);
		r.setReadBases(B("ATAA" + SequenceUtil.reverseComplement("TCTTT")));
		r.setAttribute("SA", "Contig,22,-,5M4S,0,0");
		SplitReadEvidence e = (SplitReadEvidence)asAssemblyEvidence(new MockSAMEvidenceSource(getContext(ref)), r);
		assertEquals("AAAAA", e.getHomologySequence());
		assertEquals(2, e.getHomologyAnchoredBaseCount());
	}
	public void microhomology_match_bb() {
		//          1         2         3         4
		// 1234567890123456789012345678901234567890123456789
		// AAAACCGGGCCCCAAAAAAAAAGGGGCCCGGAAAAA
		//           <<<<           <<<<<<<<<
		String contig = "AAAACCGGGCCCCAAAAAAAAAGGGGCCCGGAAAAA";
		InMemoryReferenceSequenceFile ref = new InMemoryReferenceSequenceFile(new String[] { "Contig" }, new byte[][] { B(contig) });
		SAMRecord r = new SAMRecord(new SAMFileHeader());
		r.getHeader().setSequenceDictionary(getSequenceDictionary());
		r.setReferenceIndex(0);
		r.setCigarString("9S4M");
		r.setAlignmentStart(11);
		r.setReadNegativeStrandFlag(true);
		r.setReadBases(B("TTCCGGGCCCCA"));
		r.setAttribute("SA", "Contig,36,+,4S9M,0,0");
		SplitReadEvidence e = (SplitReadEvidence)asAssemblyEvidence(new MockSAMEvidenceSource(getContext(ref)), r);
		assertEquals("CCGGGCCCC", e.getHomologySequence());
		assertEquals(3, e.getHomologyAnchoredBaseCount());
		assertEquals(11, e.getBreakendSummary().nominal);
		assertEquals(26, e.getBreakendSummary().nominal2);
		assertEquals(5, e.getBreakendSummary().start);
		assertEquals(14, e.getBreakendSummary().end);
		assertEquals(23, e.getBreakendSummary().start2);
		assertEquals(32, e.getBreakendSummary().end2);
	}
	@Test
	public void microhomology_bases_should_ignore_case() {
		SAMRecord r = withAttr("SA", "polyA,100,+,4M4S,0,0", withSequence("aaaaaaaa", Read(0, 13, "4S4M")))[0];
		SplitReadEvidence e = (SplitReadEvidence)asAssemblyEvidence(r);
		assertEquals(8, e.getHomologySequence().length());
	}
	@Test(expected=IllegalStateException.class)
	public void should_not_calculate_homology_if_reference_sequence_cannot_be_found() {
		SplitReadEvidence sr = SR(Read(0, 1, "5M5S"), Read(1, 1, "5M"));
		List<SingleReadEvidence> e = SingleReadEvidence.createEvidence(null, 0, sr.getSAMRecord());
		assertEquals(0, e.get(0).getHomologySequence().length());
	}
	@Test
	public void should_clip_breakpoints_to_contig_bounds() {
		SplitReadEvidence sr = (SplitReadEvidence)asAssemblyEvidence(withAttr("SA", "polyACGT,1,+,1X1N1X10S,0,0", Read(0, 1, "2S10M"))[0]);
		assertEquals(new BreakpointSummary(0, BWD, 1, 1, 1, 1, FWD, 2, 1, 3), sr.getBreakendSummary());
	}
	@Test
	@Category(Hg19Tests.class)
	@Ignore("outdated assembly annotations")
	public void assembly_scoring_should_be_symmetrical() throws FileNotFoundException {
		File sam = new File("src/test/resources/assembly_score_mismatch.sam");
		File ref = Hg19Tests.findHg19Reference();
		ProcessingContext pc = new ProcessingContext(getFSContext(), ref, new SynchronousReferenceLookupAdapter(new IndexedFastaSequenceFile(ref)), null, getConfig());
		SAMEvidenceSource ses = SES(pc);
		SamReader sr = SamReaderFactory.makeDefault().open(sam);
		List<SAMRecord> in = Lists.newArrayList(sr.iterator());
		List<SingleReadEvidence> sreList = in.stream()
				.flatMap(record -> SingleReadEvidence.createEvidence(ses, 0, record).stream())
				.collect(Collectors.toList());
		ImmutableListMultimap<String, SingleReadEvidence> mm = Multimaps.index(sreList, sre -> sre.getSAMRecord().getReadName());
		for (Entry<String, Collection<SingleReadEvidence>> entry : mm.asMap().entrySet()) {
			float bequal = entry.getValue().iterator().next().getBreakendQual();
			float bpqual = ((SplitReadEvidence) entry.getValue().iterator().next()).getBreakpointQual();
			for (SingleReadEvidence sre : entry.getValue()) {
				assertEquals(bequal, sre.getBreakendQual(), 0);
				assertEquals(bpqual, ((SplitReadEvidence) sre).getBreakpointQual(), 0);
			}
		}
	}
	@Test
	public void involvesPrimaryReadAlignment_should_require_not_supplementary() {
		List<SingleReadEvidence> list;
		SAMRecord r = Read(0, 1, "10M4D10M10S");
		list = SingleReadEvidence.createEvidence(SES(), 0, r);
		assertTrue(list.stream().allMatch(e -> e.involvesPrimaryReadAlignment()));
		
		r.setSupplementaryAlignmentFlag(true);
		list = SingleReadEvidence.createEvidence(SES(), 0, r);
		assertTrue(list.stream().allMatch(e -> !e.involvesPrimaryReadAlignment()));
		
		//r.setSupplementaryAlignmentFlag(false);
		//r.setNotPrimaryAlignmentFlag(true);
		//list = SingleReadEvidence.createEvidence(SES(), 0, r);
		//assertTrue(list.stream().allMatch(e -> !e.involvesPrimaryReadAlignment()));
	}
	@Test
	public void should_not_create_split_read_if_read_length_does_not_match_SA_tag() {
		SAMRecord r = withAttr("SA", "polyA,100,-,108S40M,29,1;polyA,200,-,35M113S,13,1", Read(1, 5, "36M2I10M"))[0];
		List<SingleReadEvidence> e = SingleReadEvidence.createEvidence(SES(), 0, r);
		assertEquals(0,  e.size());
	}
	@Test
	public void strandBias_should_match_read_strand() {
		assertEquals(1, SingleReadEvidence.createEvidence(SES(), 0, Read(0, 1, "50M50S")).get(0).getStrandBias(), 0);
		assertEquals(0, SingleReadEvidence.createEvidence(SES(), 0, onNegative(Read(0, 1, "50M50S"))[0]).get(0).getStrandBias(), 0);
	}
	@Test
	public void getBreakendReadOffsetInterval_should_consider_mapping_strand() {
		Range<Integer> r = SingleReadEvidence.createEvidence(SES(), 0, withSequence("NNNNN", Read(0, 1, "2M3S"))[0]).get(0).getBreakendAssemblyContigBreakpointInterval();
		// 01234
		// MMSSS
		assertEquals(2, (int)r.lowerEndpoint());
		assertEquals(2, (int)r.upperEndpoint());

		r = SingleReadEvidence.createEvidence(SES(), 0, withSequence("NNNNN", onNegative(Read(0, 1, "1M4S")))[0]).get(0).getBreakendAssemblyContigBreakpointInterval();
		// 43210
		// MSSSS
		assertEquals(4, (int)r.lowerEndpoint());
		assertEquals(4, (int)r.upperEndpoint());
	}
	@Test
	public void getBreakendReadOffsetInterval_should_untemplated_inserted_Sequence() {
		Range<Integer> r = SingleReadEvidence.createEvidence(SES(), 0, withAttr("SA", "polyA,10,+,3S2M,20,0", withSequence("NNNNN", Read(0, 1, "2M3S")))[0]).get(0).getBreakendAssemblyContigBreakpointInterval();
		// 01234
		// MMSSS
		assertEquals(2, (int)r.lowerEndpoint());
		assertEquals(3, (int)r.upperEndpoint());

		r = SingleReadEvidence.createEvidence(SES(), 0, withAttr("SA", "polyA,10,+,2M3S,20,0", withSequence("NNNNN", onNegative(Read(0, 1, "1M4S"))))[0]).get(0).getBreakendAssemblyContigBreakpointInterval();
		// 43210
		// MSSSS
		// SSSMM
		assertEquals(2, (int)r.lowerEndpoint());
		assertEquals(4, (int)r.upperEndpoint());
	}
	@Test
	public void getBreakendReadOffsetInterval_should_consider_homology() {
		Range<Integer> r = SingleReadEvidence.createEvidence(SES(), 0, withAttr("SA", "polyA,10,+,2S3M,20,0", withSequence("NAAANN", Read(0, 1, "2M4S")))[0]).get(0).getBreakendAssemblyContigBreakpointInterval();
		// 01234
		// MMSSS
		//  ***  homology interval
		assertEquals(1, (int)r.lowerEndpoint());
		assertEquals(4, (int)r.upperEndpoint());
	}
}
