package au.edu.wehi.idsv;

import static org.junit.Assert.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.function.BiConsumer;

import org.junit.Test;
import org.junit.experimental.categories.Category;
import org.junit.rules.TemporaryFolder;

import com.google.common.io.Files;

import au.edu.wehi.idsv.configuration.GridssConfiguration;
import au.edu.wehi.idsv.picard.InMemoryReferenceSequenceFile;
import au.edu.wehi.idsv.picard.SynchronousReferenceLookupAdapter;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.Header;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;


public class SplitReadEvidenceTest extends TestHelper {
	@Test
	public void should_use_SA_tag_for_split_read() {
		SAMRecord r = Read(2, 1, "3S2M5S");
		//                    rrXmm-----
		r.setReadBases(    B("ACCTAGAGGG"));
		r.setBaseQualities(B("1234567890"));
		r.setMappingQuality(40);
		r.setAttribute("SA", "polyA,100,+,2M8S,10,0");
		SplitReadEvidence e = SplitReadEvidence.create(SES(), r).get(0);
		assertEquals(new BreakpointSummary(2, BWD, 1, 0, FWD, 101), e.getBreakendSummary());
		assertEquals("ACC", S(e.getBreakendSequence()));
		assertEquals("123", S(e.getBreakendQuality()));
		assertEquals("TA", S(e.getAnchorSequence()));
		assertEquals("45", S(e.getAnchorQuality()));
		assertEquals("C", e.getUntemplatedSequence());
		assertEquals(40, e.getLocalMapq());
		assertEquals(10, e.getRemoteMapq());
	}
	@Test
	public void should_return_adjacent_alignments() {
		SAMRecord r = Read(2, 1, "3S2M5S");
		//                    AB-mm-CDDD
		r.setReadBases(    B("ACCTAGAGGG"));
		r.setBaseQualities(B("1234567890"));
		r.setMappingQuality(40);
		r.setAttribute("SA", "polyA,300,+,6S1M3S,30,0;polyA,200,+,1S1M8S,20,0;polyA,400,+,7S3M,40,0;polyA,100,+,1M9S,10,0;");
		//                                C                  B                     D                  A
		List<SplitReadEvidence> list = SplitReadEvidence.create(SES(), r);
		assertEquals(2, list.size());
		assertEquals(new BreakpointSummary(2, BWD, 1, 0, FWD, 200), list.get(0).getBreakendSummary());
		assertEquals("TA", S(list.get(0).getAnchorSequence()));
		assertEquals("ACC", S(list.get(0).getBreakendSequence()));
		assertEquals("C", list.get(0).getUntemplatedSequence());
		assertEquals(new BreakpointSummary(2, FWD, 2, 0, BWD, 300), list.get(1).getBreakendSummary());
		assertEquals("TA", S(list.get(1).getAnchorSequence()));
		assertEquals("GAGGG", S(list.get(1).getBreakendSequence()));
		assertEquals("G", list.get(1).getUntemplatedSequence());
	}
	@Test
	public void should_handle_negative_strand_remote() {
		SAMRecord r = Read(2, 1, "3S2M5S");
		//                    -BBmm-CC--
		r.setReadBases(    B("ACCNAGAGGG"));
		r.setBaseQualities(B("1234567890"));
		r.setMappingQuality(40);
		r.setAttribute("SA", "polyA,100,-,7S2M1S,0,0;polyA,200,-,2S2M6S,0,0");
		List<SplitReadEvidence> list = SplitReadEvidence.create(SES(), r);
		assertEquals(new BreakpointSummary(2, BWD, 1, 0, BWD, 100), list.get(0).getBreakendSummary());
		assertEquals("NA", S(list.get(0).getAnchorSequence()));
		assertEquals("ACC", S(list.get(0).getBreakendSequence()));
		assertEquals("", list.get(0).getUntemplatedSequence());
		assertEquals(new BreakpointSummary(2, FWD, 2, 0, FWD, 201), list.get(1).getBreakendSummary());
		assertEquals("NA", S(list.get(1).getAnchorSequence()));
		assertEquals("GAGGG", S(list.get(1).getBreakendSequence()));
		assertEquals("G", list.get(1).getUntemplatedSequence());
	}
	@Test
	public void should_handle_negative_strand_local() {
		SAMRecord r = Read(2, 1, "5S2M3S");
		//                    --CC-mmBB-
		r.setReadBases(    B("ACCTAGAGGG"));
		r.setBaseQualities(B("1234567890"));
		r.setMappingQuality(40);
		r.setReadNegativeStrandFlag(true);
		r.setAttribute("SA", "polyA,100,+,1S2M7S,0,0;polyA,200,+,6S2M2S,0,0");
		List<SplitReadEvidence> list = SplitReadEvidence.create(SES(), r);
		assertEquals(new BreakpointSummary(2, FWD, 2, 0, FWD, 101), list.get(0).getBreakendSummary());
		assertEquals("GA", S(list.get(0).getAnchorSequence()));
		assertEquals("GGG", S(list.get(0).getBreakendSequence()));
		assertEquals("", list.get(0).getUntemplatedSequence());
		assertEquals(new BreakpointSummary(2, BWD, 1, 0, BWD, 200), list.get(1).getBreakendSummary());
		assertEquals("GA", S(list.get(1).getAnchorSequence()));
		assertEquals("ACCTA", S(list.get(1).getBreakendSequence()));
		assertEquals("A", list.get(1).getUntemplatedSequence());
	}
	@Test
	public void should_consider_overlapping_alignments_as_microhomology() {
		SAMRecord r = Read(2, 1, "7M3S");
		// MMMMMMMSSS
		// SSMMMMMMMM
		//   ^---^ homology
		r.setReadBases(    B("ACCTAGAGGG"));
		r.setBaseQualities(B("1234567890"));
		r.setMappingQuality(40);
		r.setAttribute("SA", "polyA,100,+,2S8M,0,0");
		List<SplitReadEvidence> list = SplitReadEvidence.create(SES(), r);
		assertEquals(new BreakpointSummary(2, FWD, 7, 3, 7, 0, BWD, 104, 100, 104), list.get(0).getBreakendSummary());
		assertEquals("ACCTAGA", S(list.get(0).getAnchorSequence()));
		assertEquals("GGG", S(list.get(0).getBreakendSequence()));
		assertEquals("", list.get(0).getUntemplatedSequence());
	}
	@Test
	public void should_recognise_XNX_as_unanchored_breakend_interval() {
		// 1
		// 01234567
		// N--NATG
		// |--|     <- sequence expected to be placed somewhere in the N interval
		//      1
		//      0
		// 567890123456
		//     ATG  <- actually was placed
		//     SMM
		//   |--|   <- so we expect the breakpoint somewhere here
		// the soft clipping of A does not push out the expected breakend
		// position at the remote side - it could be soft clipped because
		// that is were the breakend actually is and we just haven't got
		// a good enough anchor in the contig to place it locally
		SAMRecord r = withSequence("NNATG", Read(1, 10, "1X2N1X3S"))[0];
		r.setAttribute("SA", "polyA,100,+,3S2M,0,0");
		List<SplitReadEvidence> list = SplitReadEvidence.create(SES(), r);
		assertEquals(1, list.size());
		SplitReadEvidence e = list.get(0);
		assertFalse(e.isBreakendExact());
		assertEquals(new BreakpointSummary(1, FWD, 12, 10, 13, 0, BWD, 99, 97, 100), e.getBreakendSummary());
		assertEquals("ATG", S(e.getBreakendSequence()));
		assertEquals("", S(e.getAnchorSequence()));
		assertEquals("A", e.getUntemplatedSequence());
	}
	@Test
	public void should_recogise_XNX_on_both_sides() {
		// supp alignment for should_recognise_XNX_as_unanchored_breakend_interval()
		SAMRecord r = withSequence("NNATG", Read(0, 100, "3S2M"))[0];
		r.setAttribute("SA", "polyACGT,10,+,1X2N1X3S,0,0");
		List<SplitReadEvidence> list = SplitReadEvidence.create(SES(), r);
		assertEquals(1, list.size());
		SplitReadEvidence e = list.get(0);
		assertFalse(e.isBreakendExact());
		assertEquals(new BreakpointSummary(0, BWD, 99, 97, 100, 1, FWD, 12, 10, 13), e.getBreakendSummary());
		assertEquals("A", S(e.getBreakendSequence()));
		assertEquals("TG", S(e.getAnchorSequence()));
		assertEquals("A", e.getUntemplatedSequence());
	}
	@Test
	public void should_consider_strand_for_XNX() {
		SAMRecord r = withSequence("NATG", Read(1, 10, "1X3S"))[0];
		r.setAttribute("SA", "polyA,100,-,3M1S,0,0");
		List<SplitReadEvidence> list = SplitReadEvidence.create(SES(), r);
		assertEquals(1, list.size());
		SplitReadEvidence e = list.get(0);
		assertFalse(e.isBreakendExact());
		assertEquals("ATG", S(e.getBreakendSequence()));
		assertEquals("", S(e.getAnchorSequence()));
		assertEquals(new BreakpointSummary(1, FWD, 10, 0, FWD, 102), e.getBreakendSummary());
	}
	@Test
	public void should_return_untemplated_sequence() {
		SAMRecord r = withSequence("ACTG", Read(0, 1, "1M3S"))[0];
		r.setAttribute("SA", "polyA,100,+,2S2M,10,0");
		SplitReadEvidence e = SplitReadEvidence.create(SES(), r).get(0);
		assertEquals("C", e.getUntemplatedSequence());
	}
	@Test
	public void should_return_untemplated_sequence_remote_bwd() {
		SAMRecord r = withSequence("ACTG", Read(0, 1, "1M3S"))[0];
		r.setAttribute("SA", "polyA,100,-,2M2S,10,0");
		SplitReadEvidence e = SplitReadEvidence.create(SES(), r).get(0);
		assertEquals("C", e.getUntemplatedSequence());
	}
	@Test
	public void should_return_remote_unanchored() {
		SAMRecord r = withSequence("ACN", Read(0, 1, "1M1D1M1S"))[0];
		r.setAttribute("SA", "polyA,10,+,2S1X,39,1");
		SplitReadEvidence e = SplitReadEvidence.create(SES(), r).get(0);
		assertNotNull(e);
	}
	@Test
	public void unanchored_should_handle_no_untemplated() {
		BiConsumer<String, String> go = (r1, r2) -> {
			SplitReadEvidence e1 = SplitReadEvidence.create(SES(), withAttr("SA", r1, Read(r2))[0]).get(0);
			SplitReadEvidence e2 = SplitReadEvidence.create(SES(), withAttr("SA", r2, Read(r1))[0]).get(0);
			assertEquals("", e1.getUntemplatedSequence());
			assertEquals("", e1.getUntemplatedSequence());
			assertFalse(e1.isBreakendExact());
			assertFalse(e2.isBreakendExact());
			assertEquals(3, e1.getAnchorSequence().length + e2.getAnchorSequence().length);
			assertEquals(3, e1.getBreakendSequence().length + e2.getBreakendSequence().length);
			assertEquals("", e1.getHomologySequence());
			assertEquals("", e2.getHomologySequence());
			assertEquals(0, e1.getHomologyAnchoredBaseCount());
			assertEquals(0, e2.getHomologyAnchoredBaseCount());
		};
		// fwd
		go.accept("polyA,10,+,3S1X,0,0", "polyA,10,+,3M1S,0,0");
		go.accept("polyA,10,+,3S1X,0,0", "polyA,10,-,1S3M,0,0");
		go.accept("polyA,10,+,3S2X,0,0", "polyA,10,+,3M2S,0,0");
		go.accept("polyA,10,+,3S2X,0,0", "polyA,10,-,2S3M,0,0");
		go.accept("polyA,10,+,1X3S,0,0", "polyA,10,+,1S3M,0,0");
		go.accept("polyA,10,+,1X3S,0,0", "polyA,10,-,3M1S,0,0");
		go.accept("polyA,10,+,2X3S,0,0", "polyA,10,+,2S3M,0,0");
		go.accept("polyA,10,+,2X3S,0,0", "polyA,10,-,3M2S,0,0");
	}
	@Test
	public void should_remove_split_reads_that_fully_align_to_either_side() {
		String refStr = "TAAATTGGAACACTATACCAAAACATTAACCAGCATAGCAGTATATAAGGTTAAACATTAAATAACCCCTGGCTTAACTAACTCTCCAATTGCACTTTCTATAAGTAATTGTTGTTTAGACTTTATTAATTCAGATGTTTCAGACATGTCTTATATACACAAGAGAATTTCATTTCTCTTT";
		String readStr = "AAATTGGAACACTATACCAAAACATTAACCAGCATAGCAGTATATAAGGTTAAACATTAAATAACCCCTGGCTTAACTAACTCTCCAATTGCACTTTCTATAAGTAATTGTTGTTTAGACTTTATTAATTC";
		//               1234567890123456
		//               MMMMMSSSSSSSSSS
		//                sssssMMMMMMMMMM
		InMemoryReferenceSequenceFile ref = new InMemoryReferenceSequenceFile(new String[] { "Contig" }, new byte[][] { B(refStr) });
		SAMRecord r = new SAMRecord(new SAMFileHeader());
		r.getHeader().setSequenceDictionary(ref.getSequenceDictionary());
		r.setReferenceIndex(0);
		r.setCigarString("97M34S");
		r.setAlignmentStart(1);
		r.setReadNegativeStrandFlag(false);
		r.setReadBases(B(readStr));
		r.setBaseQualities(B(refStr));
		SAMRecord r2 = SAMRecordUtil.realign(ref, r, 10, true);
		assertEquals(2, r2.getAlignmentStart());
		assertEquals("131M", r2.getCigarString());
	}
	@Test
	@Category(Hg19Tests.class)
	public void indel_mismapping_false_positive_assembly_should_not_throw_homology_error() throws IOException {
		File ref = Hg19Tests.findHg19Reference();
		IndexedFastaSequenceFile indexed = new IndexedFastaSequenceFile(ref);
		TemporaryFolder folder = new TemporaryFolder();
		folder.create();
		ProcessingContext pc = new ProcessingContext(
				new FileSystemContext(folder.getRoot(), folder.getRoot(), 500000), ref, new SynchronousReferenceLookupAdapter(indexed), new ArrayList<Header>(),
				getConfig());
		//File bam = new File(folder.getRoot(), "input.bam");
		//Files.copy(new File("src/test/resources/indel_mismapping_false_positive_assembly.sv.bam"), bam);
		SAMEvidenceSource ses = new MockSAMEvidenceSource(pc);
		List<SAMRecord> records = IntermediateFilesTest.getRecords(new File("src/test/resources/indel_mismapping_false_positive_assembly.sv.bam"));
		for (SAMRecord r : records) {
			if ("asm5".equals(r.getReadName())) {
				for (SplitReadEvidence e : SplitReadEvidence.create(ses, r)) {
					assertTrue(e.isReference());
				}
			}
		}
		folder.delete();
	}
}
