package au.edu.wehi.idsv;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.function.BiConsumer;

import htsjdk.samtools.*;
import org.junit.Assert;
import org.junit.Test;
import org.junit.experimental.categories.Category;
import org.junit.rules.TemporaryFolder;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;

import au.edu.wehi.idsv.picard.InMemoryReferenceSequenceFile;
import au.edu.wehi.idsv.picard.SynchronousReferenceLookupAdapter;
import au.edu.wehi.idsv.sam.ChimericAlignment;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.sam.SamTags;
import htsjdk.samtools.metrics.Header;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

import static org.junit.Assert.*;


public class SplitReadEvidenceTest extends TestHelper {
	@Test
	public void should_use_SA_tag_for_split_read() {
		SAMRecord r = Read(2, 1, "3S2M5S");
		//                    rrXmm-----
		r.setReadBases(B("ACCTAGAGGG"));
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
		r.setReadBases(B("ACCTAGAGGG"));
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
		r.setReadBases(B("ACCNAGAGGG"));
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
		r.setReadBases(B("ACCTAGAGGG"));
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
		// 1234567890
		//  |>   |>
		// MMMMMMMSSS
		// SSMMMMMMMM
		//   ^---^ homology
		//  10
		//   01234567
		//  <|   <|
		//
		r.setReadBases(B("ACCTAGAGGG"));
		r.setBaseQualities(B("1234567890"));
		r.setMappingQuality(40);
		r.setAttribute("SA", "polyA,100,+,2S8M,0,0");
		List<SplitReadEvidence> list = SplitReadEvidence.create(SES(), r);
		assertEquals(new BreakpointSummary(2, FWD, 7, 2, 7, 0, BWD, 105, 100, 105), list.get(0).getBreakendSummary());
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
		InMemoryReferenceSequenceFile ref = new InMemoryReferenceSequenceFile(new String[]{"Contig"}, new byte[][]{B(refStr)});
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

	@Test
	public void score_should_be_symmetrical_even_if_anchors_overlap() {
		SAMRecord primary = Read(0, 100, "6M4S");
		primary.setReadBases(B("NNNNNNNNNN"));
		primary.setBaseQualities(B("1234567890"));
		primary.setMappingQuality(40);

		SAMRecord supp = Read(0, 200, "3S7M");
		supp.setSupplementaryAlignmentFlag(true);
		supp.setReadBases(B("NNNNNNNNNN"));
		supp.setBaseQualities(B("1234567890"));
		supp.setMappingQuality(20);

		primary.setAttribute("SA", new ChimericAlignment(supp).toString());
		supp.setAttribute("SA", new ChimericAlignment(primary).toString());

		StubSAMEvidenceSource ses = new StubSAMEvidenceSource(getContext(), null, 0, 0, 100);
		SplitReadEvidence ep = SplitReadEvidence.create(ses, primary).get(0);
		SplitReadEvidence es = SplitReadEvidence.create(ses, supp).get(0);
		assertEquals(ep.getBreakpointQual(), es.getBreakpointQual(), 0);
	}

	@Test
	public void involvesPrimaryReadAlignment_should_allow_primary_local() {
		List<SingleReadEvidence> list;
		SAMRecord r = Read(0, 1, "10M10S");
		r.setAttribute("SA", new ChimericAlignment(Read(0, 10, "10S10M")).toString());

		list = SingleReadEvidence.createEvidence(SES(), 0, r);
		assertTrue(list.stream().allMatch(e -> e.involvesPrimaryReadAlignment()));
	}

	@Test
	public void involvesPrimaryReadAlignment_should_allow_primary_remote() {
		List<SingleReadEvidence> list;
		SAMRecord r = Read(0, 1, "10M10S");
		r.setSupplementaryAlignmentFlag(true);
		r.setAttribute("SA", new ChimericAlignment(Read(0, 10, "10S10M")).toString());

		list = SingleReadEvidence.createEvidence(SES(), 0, r);
		assertTrue(list.stream().allMatch(e -> e.involvesPrimaryReadAlignment()));
	}

	@Test
	public void involvesPrimaryReadAlignment_should_require_either_primary() {
		List<SingleReadEvidence> list;
		SAMRecord r = Read(0, 1, "10M10S");
		r.setSupplementaryAlignmentFlag(true);
		r.setAttribute("SA", "polyA,20,+,15S5M,39,1;polyA,10,+,10S5M5S,39,1");
		// primary alignment is 15S5M at polyA:20 which does not involve this read

		list = SingleReadEvidence.createEvidence(SES(), 0, r);
		assertTrue(list.stream().allMatch(e -> !e.involvesPrimaryReadAlignment()));
	}

	@Test
	public void isReference_should_be_true_if_either_alignment_could_be_moved_to_the_other_breakend() {
		SAMRecord r1 = withSequence("TTTAACAA", Read(0, 10, "3M5S"))[0];
		SAMRecord r2 = withSequence("TTTAACAA", Read(0, 20, "3S2M3S"))[0];
		r1.setAttribute("SA", new ChimericAlignment(r2).toString());
		r2.setAttribute("SA", new ChimericAlignment(r1).toString());

		SplitReadEvidence e1 = SplitReadEvidence.create(SES(), r1).get(0);
		SplitReadEvidence e2 = SplitReadEvidence.create(SES(), r2).get(0);

		Assert.assertTrue(e2.isReference());
		Assert.assertTrue(e1.isReference());
	}

	@Test
	public void isReference_should_require_entire_aligned_region_to_match() {
		SAMRecord r1 = withSequence("TAAATAAA", Read(0, 10, "3M5S"))[0];
		SAMRecord r2 = withSequence("TAAATAAA", Read(0, 20, "3S2M3S"))[0];
		r1.setAttribute("SA", new ChimericAlignment(r2).toString());
		r2.setAttribute("SA", new ChimericAlignment(r1).toString());

		SplitReadEvidence e1 = SplitReadEvidence.create(SES(), r1).get(0);
		SplitReadEvidence e2 = SplitReadEvidence.create(SES(), r2).get(0);

		Assert.assertFalse(e2.isReference());
		Assert.assertFalse(e1.isReference());
	}

	@Test
	public void getBreakpointQual_should_be_symmetrical() {
		SAMRecord r = withMapq(4, Read(0, 100, "1S2M3S"))[0];
		SAMRecord left = withMapq(5, Read(0, 100, "1M5S"))[0];
		SAMRecord right = withMapq(6, Read(0, 100, "3S3M"))[0];
		left.setSupplementaryAlignmentFlag(true);
		right.setSupplementaryAlignmentFlag(true);
		SAMRecordUtil.calculateTemplateTags(ImmutableList.of(r, left, right), ImmutableSet.of("SA"), false, false, false, false);

		List<SingleReadEvidence> re = SingleReadEvidence.createEvidence(SES(), 0, r);
		List<SingleReadEvidence> lefte = SingleReadEvidence.createEvidence(SES(), 0, left);
		List<SingleReadEvidence> righte = SingleReadEvidence.createEvidence(SES(), 0, right);
		Assert.assertEquals(re.get(0).getBreakendQual(), lefte.get(0).getBreakendQual(), 0);
		Assert.assertEquals(re.get(1).getBreakendQual(), righte.get(0).getBreakendQual(), 0);
	}

	@Test
	public void overaligned_breakpoints_should_be_symmetrical() {
		// 12345
		// MMMM
		//   MMM
		// overlap = 2
		SAMRecord r1 = withSequence("NNNNN", Read(0, 100, "4M1S"))[0];
		SAMRecord r2 = withSequence("NNNNN", Read(0, 200, "2S3M"))[0];
		r1.setAttribute("SA", new ChimericAlignment(r2).toString());
		r2.setAttribute("SA", new ChimericAlignment(r1).toString());
		SplitReadEvidence e1 = (SplitReadEvidence) SingleReadEvidence.createEvidence(SES(), 0, r1).get(0);
		SplitReadEvidence e2 = (SplitReadEvidence) SingleReadEvidence.createEvidence(SES(), 0, r2).get(0);
		Assert.assertEquals(2, e2.getBreakendSummary().end - e2.getBreakendSummary().start);
		Assert.assertEquals(2, e1.getBreakendSummary().end - e1.getBreakendSummary().start);
		Assert.assertEquals(e1.getBreakendSummary(), e2.getBreakendSummary().remoteBreakpoint());
		Assert.assertEquals(e2.getBreakendSummary(), e1.getBreakendSummary().remoteBreakpoint());
	}

	@Test
	public void should_handle_multiple_overlapping_fragments() {
		SAMRecord r1 = Read(0, 100, "48M52S");
		SAMRecord r2 = Read(0, 100, "32S48M20S");
		SAMRecord r3 = Read(0, 100, "64S36M");
		r1.setReadName("R");
		r2.setReadName("R");
		r3.setReadName("R");
		r1.setAttribute("SA", new ChimericAlignment(r2).toString() + ";" + new ChimericAlignment(r3).toString());
		r2.setAttribute("SA", new ChimericAlignment(r1).toString() + ";" + new ChimericAlignment(r3).toString());
		r3.setAttribute("SA", new ChimericAlignment(r1).toString() + ";" + new ChimericAlignment(r2).toString());

		List<SplitReadEvidence> e1 = SplitReadEvidence.create(SES(), r1);
		List<SplitReadEvidence> e2 = SplitReadEvidence.create(SES(), r2);
		List<SplitReadEvidence> e3 = SplitReadEvidence.create(SES(), r3);
		Assert.assertEquals(1, e1.size());
		Assert.assertEquals(2, e2.size());
		Assert.assertEquals(1, e3.size());
		Assert.assertEquals(e1.get(0).getEvidenceID(), e2.get(0).getRemoteEvidenceID());
		Assert.assertEquals(e2.get(0).getEvidenceID(), e1.get(0).getRemoteEvidenceID());
		Assert.assertEquals(e2.get(1).getEvidenceID(), e3.get(0).getRemoteEvidenceID());
		Assert.assertEquals(e3.get(0).getEvidenceID(), e2.get(1).getRemoteEvidenceID());
	}

	@Test
	public void untemplated_inserted_sequence_should_report_max_over_interval() {
		SAMRecord primary = withSequence("NNNNNN", Read(1, 200, "1M5S"))[0];
		primary.setMappingQuality(100);
		SAMRecord r = Read(1, 100, "2M4S");
		r.setMappingQuality(100);
		r.setReadNegativeStrandFlag(true);
		r.setAttribute("SA", new ChimericAlignment(primary).toString());
		r.setAttribute(SamTags.IS_ASSEMBLY, 1);
		r.setAttribute(SamTags.ASSEMBLY_DIRECTION, "f");
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_TYPE, new byte[]{0, 1, 1});
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_CATEGORY, new int[]{0, 0, 1});
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_OFFSET_START, new int[]{0, 2, 3});
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_OFFSET_END, new int[]{6, 3, 5});
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_QUAL, new float[]{10, 20, 5});
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_EVIDENCEID, "1 2 3");
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_FRAGMENTID, "1 2 3");
		AssemblyEvidenceSource aes = new MockAssemblyEvidenceSource(getContext(), ImmutableList.of(SES(0), SES(1)), new File("test.bam"));
		SplitReadEvidence e = SplitReadEvidence.create(aes, r).get(0);
		// 0 1 2 3 4 5 6
		//  M S S S    primary
		//  S S S S M M realigned
		Assert.assertEquals(10, e.getBreakpointQual(), 0);
	}

	@Test
	public void homology_should_report_min_over_interval() {
		SAMRecord primary = withSequence("NAAAAN", Read(1, 200, "1M5S"))[0];
		primary.setMappingQuality(100);
		SAMRecord r = withSequence("NAAAAN", Read(1, 300, "1S5M"))[0];
		r.setMappingQuality(100);
		r.setAttribute("SA", new ChimericAlignment(primary).toString());
		r.setAttribute(SamTags.IS_ASSEMBLY, 1);
		r.setAttribute(SamTags.ASSEMBLY_DIRECTION, "f");
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_TYPE, new byte[]{0, 1, 1});
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_CATEGORY, new int[]{0, 0, 1});
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_OFFSET_START, new int[]{3, 2, 3});
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_OFFSET_END, new int[]{4, 3, 5});
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_QUAL, new float[]{10, 20, 5});
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_EVIDENCEID, "1 2 3");
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_FRAGMENTID, "1 2 3");
		AssemblyEvidenceSource aes = new MockAssemblyEvidenceSource(getContext(), ImmutableList.of(SES(0), SES(1)), new File("test.bam"));
		SplitReadEvidence e = SplitReadEvidence.create(aes, r).get(0);
		Assert.assertEquals(0, e.getBreakpointQual(), 0);
	}

	@Test
	public void FWD_unanchored_assembly_should_pro_rata_support_to_start_of_contig() {
		SAMRecord r = withSequence("NNACTG", Read(1, 100, "1X10N1X4S"))[0];
		SAMRecord realigned = Read(2, 100, "2S4M");
		r.setMappingQuality(60);
		realigned.setMappingQuality(60);
		r.setAttribute("SA", new ChimericAlignment(realigned).toString());
		realigned.setAttribute("SA", new ChimericAlignment(r).toString());
		r.setAttribute(SamTags.IS_ASSEMBLY, 1);
		r.setAttribute(SamTags.ASSEMBLY_DIRECTION, 'f');
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_TYPE, new byte[]{1});
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_CATEGORY, new int[]{0});
		// 0 1 2 3 4 5 6
		//  N N A C T G
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_OFFSET_START, new int[]{2});
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_OFFSET_END, new int[]{5});
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_QUAL, new float[]{10});
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_EVIDENCEID, "e");
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_FRAGMENTID, "e");
		r.setAttribute(SamTags.UNANCHORED, 1);
		copyAssemblyAttributes(r, realigned);
		AssemblyEvidenceSource aes = new MockAssemblyEvidenceSource(getContext(), ImmutableList.of(SES(0), SES(1)), new File("test.bam"));
		SplitReadEvidence local = SplitReadEvidence.create(aes, r).get(0);
		SplitReadEvidence remote = SplitReadEvidence.create(aes, realigned).get(0);
		Assert.assertEquals(10, local.getBreakpointQual(), 0);
		Assert.assertEquals(10, remote.getBreakpointQual(), 0);
	}

	@Test
	public void BWD_unanchored_assembly_should_pro_rata_support_to_start_of_contig() {
		SAMRecord r = withSequence("ACTGNN", Read(1, 100, "4S1X10N1X"))[0];
		SAMRecord realigned = Read(2, 100, "4M2S");
		r.setMappingQuality(60);
		realigned.setMappingQuality(60);
		r.setAttribute("SA", new ChimericAlignment(realigned).toString());
		realigned.setAttribute("SA", new ChimericAlignment(r).toString());
		r.setAttribute(SamTags.IS_ASSEMBLY, 1);
		r.setAttribute(SamTags.ASSEMBLY_DIRECTION, 'b');
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_TYPE, new byte[]{1});
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_CATEGORY, new int[]{0});
		// 0 1 2 3 4 5 6
		//  A C T G N N
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_OFFSET_START, new int[]{1});
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_OFFSET_END, new int[]{4});
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_QUAL, new float[]{10});
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_EVIDENCEID, "e");
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_FRAGMENTID, "e");
		r.setAttribute(SamTags.UNANCHORED, 1);
		copyAssemblyAttributes(r, realigned);
		AssemblyEvidenceSource aes = new MockAssemblyEvidenceSource(getContext(), ImmutableList.of(SES(0), SES(1)), new File("test.bam"));
		SplitReadEvidence local = SplitReadEvidence.create(aes, r).get(0);
		SplitReadEvidence remote = SplitReadEvidence.create(aes, realigned).get(0);
		Assert.assertEquals(10, local.getBreakpointQual(), 0);
		Assert.assertEquals(10, remote.getBreakpointQual(), 0);
	}

	@Test
	public void assembly_prorata_with_inserted_sequence_should_be_symmetrical() {

		SAMRecord r = Read(2, 300, "10S5M");
		SAMRecord r2 = Read(2, 200, "5M10S");
		r.setMappingQuality(100);
		r2.setMappingQuality(100);
		r.setAttribute("SA", new ChimericAlignment(r2).toString());
		r2.setAttribute("SA", new ChimericAlignment(r).toString());
		r.setAttribute(SamTags.IS_ASSEMBLY, 1);
		r.setAttribute(SamTags.ASSEMBLY_DIRECTION, "f");
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_TYPE, new byte[]{0, 0, 0});
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_CATEGORY, new int[]{0, 0, 0});
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_OFFSET_START, new int[]{2, 6, 0});
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_OFFSET_END, new int[]{7, 15, 15});
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_QUAL, new float[]{1, 2, 3});
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_EVIDENCEID, "1 2 3");
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_FRAGMENTID, "1 2 3");
		copyAssemblyAttributes(r, r2);
		AssemblyEvidenceSource aes = new MockAssemblyEvidenceSource(getContext(), ImmutableList.of(SES(0), SES(1)), new File("test.bam"));
		SplitReadEvidence elocal = SplitReadEvidence.create(aes, r2).get(0);
		SplitReadEvidence eremote = SplitReadEvidence.create(aes, r).get(0);
		Assert.assertEquals(elocal.getBreakpointQual(), eremote.getBreakpointQual(), 0);
	}

	public void copyAssemblyAttributes(SAMRecord source, SAMRecord dest) {
		dest.setAttribute(SamTags.IS_ASSEMBLY, source.getAttribute(SamTags.IS_ASSEMBLY));
		dest.setAttribute(SamTags.ASSEMBLY_EVIDENCE_TYPE, source.getAttribute(SamTags.ASSEMBLY_EVIDENCE_TYPE));
		dest.setAttribute(SamTags.ASSEMBLY_EVIDENCE_CATEGORY, source.getAttribute(SamTags.ASSEMBLY_EVIDENCE_CATEGORY));
		dest.setAttribute(SamTags.ASSEMBLY_EVIDENCE_OFFSET_START, source.getAttribute(SamTags.ASSEMBLY_EVIDENCE_OFFSET_START));
		dest.setAttribute(SamTags.ASSEMBLY_EVIDENCE_OFFSET_END, source.getAttribute(SamTags.ASSEMBLY_EVIDENCE_OFFSET_END));
		dest.setAttribute(SamTags.ASSEMBLY_EVIDENCE_QUAL, source.getAttribute(SamTags.ASSEMBLY_EVIDENCE_QUAL));
		dest.setAttribute(SamTags.ASSEMBLY_EVIDENCE_EVIDENCEID, source.getAttribute(SamTags.ASSEMBLY_EVIDENCE_EVIDENCEID));
		dest.setAttribute(SamTags.ASSEMBLY_EVIDENCE_FRAGMENTID, source.getAttribute(SamTags.ASSEMBLY_EVIDENCE_FRAGMENTID));
		dest.setAttribute(SamTags.UNANCHORED, source.getAttribute(SamTags.UNANCHORED));
	}

	@Test
	public void assembly_prorata_should_handle_truncated_split_read_assemblies() {
		String primary = "asm44-626172\t0\tpolyA\t200\t70\t169S2X\t*\t0\t0\tTTGACTTGTATGTCCAATAAAACACAAATTAAAACTTGAGAGAAACTTGTGAGATGTCTGCAGCACTTTAGGCTGGTTGGGACTATGAAGAAGAAATGTGTCAATTTTCTAAGCCTGTAAAGCAGTATATTCGGTGAACCCTAAAATGAATTTCACCCTTATGAAATGGNN\tBBCDEEEEIMRX]____________________________________________________________________________________________________________________________________________________________!!\tSA:Z:polyA,100,+,149M5D20M2S,60,6\tNM:i:2\taa:i:1\tua:i:1\tsb:f:0.486486\tec:B:i,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,0,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,0,1,1,1,1\tad:A:b\toe:B:i,166,167,167,167,167,167,167,167,168,167,167,167,167,167,167,164,160,167,167,167,167,167,167,167,163,167,167,167,167,167,167,168,167,167,167,164,167,165,158,167,150,167,167,167,167,167,167,167,167,167,167,167,167,122,167,161,167,167,167,167,167,167,167,167,167,167,167,167,167,167,167,167,167,167\tef:Z:ST-E00287:121:HCT3FCCXY:5:1101:30259:3630 ST-E00287:121:HCT3FCCXY:1:1213:11424:28980 ST-E00287:121:HCT3FCCXY:2:1206:2514:12842 ST-E00287:121:HCT3FCCXY:8:1110:2940:41743 ST-E00287:121:HCT3FCCXY:5:1203:3884:15883 ST-E00287:121:HCT3FCCXY:7:2201:17391:57424 ST-E00287:121:HCT3FCCXY:1:1107:7192:54910 ST-E00287:121:HCT3FCCXY:5:2222:5365:34025 ST-E00287:121:HCT3FCCXY:4:2101:11170:72209 ST-E00287:121:HCT3FCCXY:8:1217:15432:39792 ST-E00287:121:HCT3FCCXY:2:1107:5091:66285 ST-E00287:121:HCT3FCCXY:5:2108:1397:41849 ST-E00287:121:HCT3FCCXY:3:2212:8694:60536 ST-E00287:121:HCT3FCCXY:4:2217:9485:14195 ST-E00287:121:HCT3FCCXY:7:1124:7852:36821 ST-E00287:121:HCT3FCCXY:3:2207:12966:14037 ST-E00287:121:HCT3FCCXY:1:1123:8379:70926 ST-E00287:121:HCT3FCCXY:8:1107:14336:60079 ST-E00287:121:HCT3FCCXY:4:2119:24911:62189 ST-E00287:121:HCT3FCCXY:3:1124:3488:33867 ST-E00287:121:HCT3FCCXY:4:2214:6837:38333 ST-E00287:121:HCT3FCCXY:4:1117:23236:45822 ST-E00287:121:HCT3FCCXY:1:2212:18710:17061 ST-E00287:121:HCT3FCCXY:4:1216:2189:50779 ST-E00287:121:HCT3FCCXY:2:1107:13961:15971 ST-E00287:121:HCT3FCCXY:1:2214:16386:30650 ST-E00287:121:HCT3FCCXY:7:2119:5883:70926 ST-E00287:121:HCT3FCCXY:3:1213:2077:11101 ST-E00287:121:HCT3FCCXY:6:2209:17574:42657 ST-E00287:121:HCT3FCCXY:3:2118:17350:51799 ST-E00287:121:HCT3FCCXY:5:1116:6563:25587 ST-E00287:121:HCT3FCCXY:5:2205:22485:64492 ST-E00287:121:HCT3FCCXY:3:2109:23348:70943 ST-E00287:121:HCT3FCCXY:6:2203:11667:69344 ST-E00287:121:HCT3FCCXY:2:2124:10592:32320 ST-E00287:121:HCT3FCCXY:4:2109:15808:25394 ST-E00287:121:HCT3FCCXY:8:2122:20669:56774 ST-E00287:121:HCT3FCCXY:2:1105:27722:41181 ST-E00287:121:HCT3FCCXY:4:1219:8643:32320 ST-E00287:121:HCT3FCCXY:4:1105:14154:56352 ST-E00287:121:HCT3FCCXY:2:2204:3711:58620 ST-E00287:121:HCT3FCCXY:6:1117:3457:48089 ST-E00287:121:HCT3FCCXY:1:2219:8552:67006 ST-E00287:121:HCT3FCCXY:4:2210:8420:62628 ST-E00287:121:HCT3FCCXY:3:1206:10947:64052 ST-E00287:121:HCT3FCCXY:8:1105:14834:63859 ST-E00287:121:HCT3FCCXY:1:1112:32725:40671 ST-E00287:121:HCT3FCCXY:3:2102:12012:31125 ST-E00287:121:HCT3FCCXY:7:2119:5873:70838 ST-E00287:121:HCT3FCCXY:6:1115:5294:16428 ST-E00287:121:HCT3FCCXY:5:1107:19228:46156 ST-E00287:121:HCT3FCCXY:5:2124:12185:9519 ST-E00287:121:HCT3FCCXY:7:1104:2534:17061 ST-E00287:121:HCT3FCCXY:8:1211:9171:55455 ST-E00287:121:HCT3FCCXY:5:1115:9323:13668 ST-E00287:121:HCT3FCCXY:2:2115:22719:38561 ST-E00287:121:HCT3FCCXY:8:1205:11099:14178 ST-E00287:121:HCT3FCCXY:3:1224:12327:8323 ST-E00287:121:HCT3FCCXY:1:1118:31040:61521 ST-E00287:121:HCT3FCCXY:4:1115:15432:21755 ST-E00287:121:HCT3FCCXY:7:2224:17350:47650 ST-E00287:121:HCT3FCCXY:1:2201:9323:68131 ST-E00287:121:HCT3FCCXY:4:1214:19116:39352 ST-E00287:121:HCT3FCCXY:3:2223:26585:30808 ST-E00287:121:HCT3FCCXY:5:1212:10490:22335 ST-E00287:121:HCT3FCCXY:3:1212:10541:32584 ST-E00287:121:HCT3FCCXY:3:2122:14499:47281 ST-E00287:121:HCT3FCCXY:4:1118:23419:62698 ST-E00287:121:HCT3FCCXY:4:2111:8126:33111 ST-E00287:121:HCT3FCCXY:3:1219:18325:2645 ST-E00287:121:HCT3FCCXY:1:1210:18771:54682 ST-E00287:121:HCT3FCCXY:6:1208:24332:62945 ST-E00287:121:HCT3FCCXY:7:2109:22069:20137 ST-E00287:121:HCT3FCCXY:1:1103:19735:13176\tdl:i:151\teq:B:f,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.2817,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.2817,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.2817,24.297,24.297,24.297,24.297,24.297,24.2817,24.2817,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.2817,24.297,24.297,24.297,24.297,24.2817,24.297,24.297,24.297,24.297\tos:B:i,17,39,65,91,93,125,107,23,19,113,46,93,29,70,31,116,11,24,73,75,124,30,84,25,14,48,114,115,61,38,57,19,124,107,124,15,122,16,9,31,1,109,57,85,68,47,58,71,114,28,83,91,37,10,82,12,47,107,26,125,24,48,68,95,37,83,118,82,122,22,87,70,60,44\tat:B:i,0,0\tet:B:c,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0\tez:Z:1z0bkUrw63V3fjiU9swPD5hy215JhZj5 27mACrgnqJn9k2C--LpjKJUSS7LfD3kE 2P0UfQQ_UP5OQNe63vyQtg5DSDHFe3lc 3PQnGGci1M5W7sOAq6JZXrxp49xDafo9 3uRjr-vJ2XIsNQVjHT233B4qlcyUdb6k 4JUHuBYVhAWFSWWVzbkf3DBbenaEvNxL 55KS7FS28nHKMGPv38YChcWZSoR-aGdi 5GqIzRxaNRgYnaLZsRq6fkqdDFdPd-WT 62fZd_SiLF3m-6DZLMp2cWUi99F01bSB 66GkWqGiNgQJqVvLZx_tfsVkd7DRy8FT 6xpGALGBvLYtMPcBjerkT4xT2Xu1bQnG 7pNpq4MgF3ev7njLJ603EdIfyxPzSpBd 8FxR4XVjmsVPnsBc1f7lFhqtJqY7yz-I 95Kztytoj4w6Bbua9-N2KFuY44YWCnor 999upeI4xWido1lm-p9H-RxzgMhpd0rh 9dlkAX14h7fn8z-AwQzmGFMmDyLRnSXJ BcqVIH2X8y6sk3h-QFFSMTni66NbGRod CYR0ZLKX9klQf8-cZm8HbUyiM6V6mnuS C_f2JbNeMHCnW8TpPFV9k3HSqT19o8ZR ClUzWC_JlBReo7EVU8zFMhrc5dV39iK9 Crldjh0llyhyjwK80bHu7W9Vr-GHsNhe FGS217q7TbZ4iY9t88iRcxuK27RhKVsr FjvVy0C1cRJ0glULzGWltJNAx55G0Yd9 FuylHaUV4SMHbTAm4HX-pGX5nMAiLX2g GA2IORH7hS-H0Fw08WtW6c5BwCkkal-e GhO3Nu-JGOQ_TxlbYqzrO_D6oe0ZpA7S JMezAdFfdK1r_Dvg4JXKSKvmlXJ2gBNz JNc0okUl0U0R5IYmxMVH1723oO2L7LsP KKl6iaC1nIRufK43nIFvYsO-3WOQXBfb KV4ZqtIMeeYQJ-3JpFZiEK7Thd7HgWcF KgyygNyGcUfQrMQvPPp5OcTWs0EVoSFZ LdM16vDCkKm32h8VcpBJnxW0gjL7uRsM LgEnMJU9we9K43f1p0Vzwhq3aGM2dJ2T NE3ajMf8nCwKIFUFIJs3JcQvcLir4m57 O6as4u54uPKOnVZw-u2SgJzUGdQ3ZUIS O8VMwwMLLRccABtsyYIGLgVOTeuWsa_d OFNmbKXUJAT16793ofOuWXdN_U_nMk6Q Obo6ZqjwzQW_hICodnCfXfNDTgak-pzq OxJFmTIL4Qssz1Huw3N4JBQYxmSSGnO1 P4LXy7kKXS17gmbBZimN10Fs12kD0KzI POrB2b1cz-uOm50XNA8MzKq-UcycfDY6 PdCmwwkOBrict3hAe9HGmcGGfyM7DDFm ULKrTLEAmiM6i5RevIktGaFyrAwhdGjV VbXIdWW_Onq0apYhynyqQ5PKNtSV8yc7 WaScPY0TGlx5WcMMwxVv0zdK6uIybBGf YqpmqFh0MGckBzk19J8QNLcN0Z57iBGh Z_rcKmc3vSMg4zZxDVkQsAzqXCMQ47Oj _vvBlG1feQouqY63554YvXuuBOAOedra a0Ya6hgRJ4UTRkk1DlUlFngWsFCx15GN cVgsVl4L0lDpHVLO6MasGSfx73S6oUd1 c_PwIuDe6Z7ZsUQc8H_Qmgzk_n-0-3Qd fQFHMqDRmLGvF4BOBoJBEo6ZN3LbLRZR h03pCdDZTQGZbvt3HbNiCOgOEo7sYtgJ h61M-LHz_9GWq2GCnw2zCDMIdOcnDXri i0Ivxte0Q8L3lB6WfroJe3WyR_I-hNnz ib0WWc6eUzdhcTQMCOJTlwMOAdNX3Xj_ ip154wrXXSn27Qm895XNkDX5AEgko3NN j4UBlNCAm4hTRSfXEYVN4jZ337M8BXNd lFjexdaOMWqzcWYIh3sJsiNDyVsRCZFe lOzjz4vRxY3_FZHG2kGaKKWC4d8TUIQ4 mGWe3EC1LyJB44aruBmBB1hprrm4reym mZJIbxLz8Dhbz5ByF-FAbPiIU5B_dfgC mozJJ6ZiDR3PUDF5YjVv1AWLwRUe2ybp mrbkhtbS1PYEU5fxoH2_tvHMVgE5Sr3G n0XognU0vW7bT_XS8hwJ5aU-y2TtGaDJ qLLyGnfNDviPa-vIxXqVxSYz6apL3ORW ssps_TvWFICPc14XyLlqLeDMhEjoJaoD ta-OGsr0F55Xtm-SJ9Uwl2_xQqRGyp1l vDFWpE9M7yfHSuIu6lGJpuuYjtuv4zgC vYpV_H8QzoSCtdpTf0xVAy0g1CZIdb_3 wBEz-ouSzeXL6VjYCxd9EEBNXAKV2Wmj wf-KITG5zs8F9igxU5kcW9e-0JMeZYyO ydeOuqNKJCUhGjAbjmt4tePKOwkzsy13 zU6-mWnFMhB_XubGtq9ikTHeY2xIHy1u\n";
		String supp = "asm44-626172\t2048\tpolyA\t100\t60\t149M5D20M2S\t*\t0\t0\tTTGACTTGTATGTCCAATAAAACACAAATTAAAACTTGAGAGAAACTTGTGAGATGTCTGCAGCACTTTAGGCTGGTTGGGACTATGAAGAAGAAATGTGTCAATTTTCTAAGCCTGTAAAGCAGTATATTCGGTGAACCCTAAAATGAATTTCACCCTTATGAAATGGNN\tBBCDEEEEIMRX]____________________________________________________________________________________________________________________________________________________________!!\tZSA::polyA,200,+,169S2X,70,2\tNM:i:6\tAS:i:157\tXS:i:22\taa:i:1\tua:i:1\tsb:f:0.486486\tec:B:i,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,0,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,0,1,1,1,1\tad:A:b\toe:B:i,166,167,167,167,167,167,167,167,168,167,167,167,167,167,167,164,160,167,167,167,167,167,167,167,163,167,167,167,167,167,167,168,167,167,167,164,167,165,158,167,150,167,167,167,167,167,167,167,167,167,167,167,167,122,167,161,167,167,167,167,167,167,167,167,167,167,167,167,167,167,167,167,167,167\tef:Z:ST-E00287:121:HCT3FCCXY:5:1101:30259:3630 ST-E00287:121:HCT3FCCXY:1:1213:11424:28980 ST-E00287:121:HCT3FCCXY:2:1206:2514:12842 ST-E00287:121:HCT3FCCXY:8:1110:2940:41743 ST-E00287:121:HCT3FCCXY:5:1203:3884:15883 ST-E00287:121:HCT3FCCXY:7:2201:17391:57424 ST-E00287:121:HCT3FCCXY:1:1107:7192:54910 ST-E00287:121:HCT3FCCXY:5:2222:5365:34025 ST-E00287:121:HCT3FCCXY:4:2101:11170:72209 ST-E00287:121:HCT3FCCXY:8:1217:15432:39792 ST-E00287:121:HCT3FCCXY:2:1107:5091:66285 ST-E00287:121:HCT3FCCXY:5:2108:1397:41849 ST-E00287:121:HCT3FCCXY:3:2212:8694:60536 ST-E00287:121:HCT3FCCXY:4:2217:9485:14195 ST-E00287:121:HCT3FCCXY:7:1124:7852:36821 ST-E00287:121:HCT3FCCXY:3:2207:12966:14037 ST-E00287:121:HCT3FCCXY:1:1123:8379:70926 ST-E00287:121:HCT3FCCXY:8:1107:14336:60079 ST-E00287:121:HCT3FCCXY:4:2119:24911:62189 ST-E00287:121:HCT3FCCXY:3:1124:3488:33867 ST-E00287:121:HCT3FCCXY:4:2214:6837:38333 ST-E00287:121:HCT3FCCXY:4:1117:23236:45822 ST-E00287:121:HCT3FCCXY:1:2212:18710:17061 ST-E00287:121:HCT3FCCXY:4:1216:2189:50779 ST-E00287:121:HCT3FCCXY:2:1107:13961:15971 ST-E00287:121:HCT3FCCXY:1:2214:16386:30650 ST-E00287:121:HCT3FCCXY:7:2119:5883:70926 ST-E00287:121:HCT3FCCXY:3:1213:2077:11101 ST-E00287:121:HCT3FCCXY:6:2209:17574:42657 ST-E00287:121:HCT3FCCXY:3:2118:17350:51799 ST-E00287:121:HCT3FCCXY:5:1116:6563:25587 ST-E00287:121:HCT3FCCXY:5:2205:22485:64492 ST-E00287:121:HCT3FCCXY:3:2109:23348:70943 ST-E00287:121:HCT3FCCXY:6:2203:11667:69344 ST-E00287:121:HCT3FCCXY:2:2124:10592:32320 ST-E00287:121:HCT3FCCXY:4:2109:15808:25394 ST-E00287:121:HCT3FCCXY:8:2122:20669:56774 ST-E00287:121:HCT3FCCXY:2:1105:27722:41181 ST-E00287:121:HCT3FCCXY:4:1219:8643:32320 ST-E00287:121:HCT3FCCXY:4:1105:14154:56352 ST-E00287:121:HCT3FCCXY:2:2204:3711:58620 ST-E00287:121:HCT3FCCXY:6:1117:3457:48089 ST-E00287:121:HCT3FCCXY:1:2219:8552:67006 ST-E00287:121:HCT3FCCXY:4:2210:8420:62628 ST-E00287:121:HCT3FCCXY:3:1206:10947:64052 ST-E00287:121:HCT3FCCXY:8:1105:14834:63859 ST-E00287:121:HCT3FCCXY:1:1112:32725:40671 ST-E00287:121:HCT3FCCXY:3:2102:12012:31125 ST-E00287:121:HCT3FCCXY:7:2119:5873:70838 ST-E00287:121:HCT3FCCXY:6:1115:5294:16428 ST-E00287:121:HCT3FCCXY:5:1107:19228:46156 ST-E00287:121:HCT3FCCXY:5:2124:12185:9519 ST-E00287:121:HCT3FCCXY:7:1104:2534:17061 ST-E00287:121:HCT3FCCXY:8:1211:9171:55455 ST-E00287:121:HCT3FCCXY:5:1115:9323:13668 ST-E00287:121:HCT3FCCXY:2:2115:22719:38561 ST-E00287:121:HCT3FCCXY:8:1205:11099:14178 ST-E00287:121:HCT3FCCXY:3:1224:12327:8323 ST-E00287:121:HCT3FCCXY:1:1118:31040:61521 ST-E00287:121:HCT3FCCXY:4:1115:15432:21755 ST-E00287:121:HCT3FCCXY:7:2224:17350:47650 ST-E00287:121:HCT3FCCXY:1:2201:9323:68131 ST-E00287:121:HCT3FCCXY:4:1214:19116:39352 ST-E00287:121:HCT3FCCXY:3:2223:26585:30808 ST-E00287:121:HCT3FCCXY:5:1212:10490:22335 ST-E00287:121:HCT3FCCXY:3:1212:10541:32584 ST-E00287:121:HCT3FCCXY:3:2122:14499:47281 ST-E00287:121:HCT3FCCXY:4:1118:23419:62698 ST-E00287:121:HCT3FCCXY:4:2111:8126:33111 ST-E00287:121:HCT3FCCXY:3:1219:18325:2645 ST-E00287:121:HCT3FCCXY:1:1210:18771:54682 ST-E00287:121:HCT3FCCXY:6:1208:24332:62945 ST-E00287:121:HCT3FCCXY:7:2109:22069:20137 ST-E00287:121:HCT3FCCXY:1:1103:19735:13176\tdl:i:151\teq:B:f,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.2817,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.2817,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.2817,24.297,24.297,24.297,24.297,24.297,24.2817,24.2817,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.297,24.2817,24.297,24.297,24.297,24.297,24.2817,24.297,24.297,24.297,24.297\tos:B:i,17,39,65,91,93,125,107,23,19,113,46,93,29,70,31,116,11,24,73,75,124,30,84,25,14,48,114,115,61,38,57,19,124,107,124,15,122,16,9,31,1,109,57,85,68,47,58,71,114,28,83,91,37,10,82,12,47,107,26,125,24,48,68,95,37,83,118,82,122,22,87,70,60,44\tat:B:i,0,0\tet:B:c,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0\tez:Z:1z0bkUrw63V3fjiU9swPD5hy215JhZj5 27mACrgnqJn9k2C--LpjKJUSS7LfD3kE 2P0UfQQ_UP5OQNe63vyQtg5DSDHFe3lc 3PQnGGci1M5W7sOAq6JZXrxp49xDafo9 3uRjr-vJ2XIsNQVjHT233B4qlcyUdb6k 4JUHuBYVhAWFSWWVzbkf3DBbenaEvNxL 55KS7FS28nHKMGPv38YChcWZSoR-aGdi 5GqIzRxaNRgYnaLZsRq6fkqdDFdPd-WT 62fZd_SiLF3m-6DZLMp2cWUi99F01bSB 66GkWqGiNgQJqVvLZx_tfsVkd7DRy8FT 6xpGALGBvLYtMPcBjerkT4xT2Xu1bQnG 7pNpq4MgF3ev7njLJ603EdIfyxPzSpBd 8FxR4XVjmsVPnsBc1f7lFhqtJqY7yz-I 95Kztytoj4w6Bbua9-N2KFuY44YWCnor 999upeI4xWido1lm-p9H-RxzgMhpd0rh 9dlkAX14h7fn8z-AwQzmGFMmDyLRnSXJ BcqVIH2X8y6sk3h-QFFSMTni66NbGRod CYR0ZLKX9klQf8-cZm8HbUyiM6V6mnuS C_f2JbNeMHCnW8TpPFV9k3HSqT19o8ZR ClUzWC_JlBReo7EVU8zFMhrc5dV39iK9 Crldjh0llyhyjwK80bHu7W9Vr-GHsNhe FGS217q7TbZ4iY9t88iRcxuK27RhKVsr FjvVy0C1cRJ0glULzGWltJNAx55G0Yd9 FuylHaUV4SMHbTAm4HX-pGX5nMAiLX2g GA2IORH7hS-H0Fw08WtW6c5BwCkkal-e GhO3Nu-JGOQ_TxlbYqzrO_D6oe0ZpA7S JMezAdFfdK1r_Dvg4JXKSKvmlXJ2gBNz JNc0okUl0U0R5IYmxMVH1723oO2L7LsP KKl6iaC1nIRufK43nIFvYsO-3WOQXBfb KV4ZqtIMeeYQJ-3JpFZiEK7Thd7HgWcF KgyygNyGcUfQrMQvPPp5OcTWs0EVoSFZ LdM16vDCkKm32h8VcpBJnxW0gjL7uRsM LgEnMJU9we9K43f1p0Vzwhq3aGM2dJ2T NE3ajMf8nCwKIFUFIJs3JcQvcLir4m57 O6as4u54uPKOnVZw-u2SgJzUGdQ3ZUIS O8VMwwMLLRccABtsyYIGLgVOTeuWsa_d OFNmbKXUJAT16793ofOuWXdN_U_nMk6Q Obo6ZqjwzQW_hICodnCfXfNDTgak-pzq OxJFmTIL4Qssz1Huw3N4JBQYxmSSGnO1 P4LXy7kKXS17gmbBZimN10Fs12kD0KzI POrB2b1cz-uOm50XNA8MzKq-UcycfDY6 PdCmwwkOBrict3hAe9HGmcGGfyM7DDFm ULKrTLEAmiM6i5RevIktGaFyrAwhdGjV VbXIdWW_Onq0apYhynyqQ5PKNtSV8yc7 WaScPY0TGlx5WcMMwxVv0zdK6uIybBGf YqpmqFh0MGckBzk19J8QNLcN0Z57iBGh Z_rcKmc3vSMg4zZxDVkQsAzqXCMQ47Oj _vvBlG1feQouqY63554YvXuuBOAOedra a0Ya6hgRJ4UTRkk1DlUlFngWsFCx15GN cVgsVl4L0lDpHVLO6MasGSfx73S6oUd1 c_PwIuDe6Z7ZsUQc8H_Qmgzk_n-0-3Qd fQFHMqDRmLGvF4BOBoJBEo6ZN3LbLRZR h03pCdDZTQGZbvt3HbNiCOgOEo7sYtgJ h61M-LHz_9GWq2GCnw2zCDMIdOcnDXri i0Ivxte0Q8L3lB6WfroJe3WyR_I-hNnz ib0WWc6eUzdhcTQMCOJTlwMOAdNX3Xj_ ip154wrXXSn27Qm895XNkDX5AEgko3NN j4UBlNCAm4hTRSfXEYVN4jZ337M8BXNd lFjexdaOMWqzcWYIh3sJsiNDyVsRCZFe lOzjz4vRxY3_FZHG2kGaKKWC4d8TUIQ4 mGWe3EC1LyJB44aruBmBB1hprrm4reym mZJIbxLz8Dhbz5ByF-FAbPiIU5B_dfgC mozJJ6ZiDR3PUDF5YjVv1AWLwRUe2ybp mrbkhtbS1PYEU5fxoH2_tvHMVgE5Sr3G n0XognU0vW7bT_XS8hwJ5aU-y2TtGaDJ qLLyGnfNDviPa-vIxXqVxSYz6apL3ORW ssps_TvWFICPc14XyLlqLeDMhEjoJaoD ta-OGsr0F55Xtm-SJ9Uwl2_xQqRGyp1l vDFWpE9M7yfHSuIu6lGJpuuYjtuv4zgC vYpV_H8QzoSCtdpTf0xVAy0g1CZIdb_3 wBEz-ouSzeXL6VjYCxd9EEBNXAKV2Wmj wf-KITG5zs8F9igxU5kcW9e-0JMeZYyO ydeOuqNKJCUhGjAbjmt4tePKOwkzsy13 zU6-mWnFMhB_XubGtq9ikTHeY2xIHy1u\n";
		final SAMLineParser parser = new SAMLineParser(new DefaultSAMRecordFactory(), ValidationStringency.DEFAULT_STRINGENCY, getHeader(), null, null);
		SAMRecord r = parser.parseLine(primary);
		AssemblyAttributes aa = new AssemblyAttributes(r);
		List<SplitReadEvidence> e = SplitReadEvidence.create(SES(), r);
		assertEquals(1, e.size());
		double qual = e.get(0).getBreakpointQual();
		assertNotEquals(0, qual, 0);
	}
    @Test
    public void assembly_prorata_should_handle_soft_clip_and_indel() {
        String primary = "asm0-438196	0	polyA	10	60	107S132M5I146M	*	0	0	TAATTGTTTGACTTTTTAGTGCTACAGTTCAAGAGTTCTTTACTTGTTATAGGTACTAGTCATTTGCCAGATGTGGGTTTTGTAAATGTTTTCTCCCATTCTGTAGCTTCTACAGGCTGTGGAGGAAGTATGGCCATGAGCATCTTCTTGGGGAGGCATCAGGGAGCTTTTACTCAGGGCGGAGGCGAAGCAGGAGCTAGAGAGTGAGGGGGAGGTGCCACACACTTTCAAACAAGCAGATCTCATAACTCTCTGTTATCAGGACAGGACCAAGCGGATGGTGCTAACCTCTTCATGAGAAATCCGCCCCCATGATCCAGTCACCTCCTACCAGGCCCCTCCTCCAACACTGAGGATTACATTTCATTATGAGATGTGGGCAGGGACACA	BBCDEEEFFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFFFFFEEEDDDCCDEFGHILPV_______________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________ZE8...	OA:Z:1,3724506,+,244S146M,60,1	SA:Z:random,100,-,283S107M,60,0	NM:i:1	aa:i:1	sb:f:0.632479	ec:B:i,1,1,1,1,1,1,1,1,0,1,0,0,1,1,1,1,0,1,1,0,1,1,1,1,0,1,1,1,1,1,1,0,1,1,0,1,1,0,0,1,1,0,1,0,1,0,1,1,1,0,0,0,0,1,0,0,0,1,1,0,1,0,1,1,1,0,0,0,0,1,1,1,0,0,1,1,1,0,1,1,0,1,1,1,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,0,1,1,1,0,1,1,1,0,1,1,1,1,0,0,1,0	ad:A:b	oe:B:i,341,243,362,271,326,335,317,310,354,299,332,355,376,335,373,345,342,383,286,341,363,296,340,356,301,265,281,386,385,283,317,305,376,348,384,285,382,345,366,262,351,362,338,355,386,377,271,321,262,371,278,304,376,261,323,253,283,300,358,377,292,308,364,384,293,369,309,293,258,269,296,303,310,355,313,365,351,370,295,335,379,226,368,369,310,370,315,323,263,318,322,263,376,300,379,150,370,372,263,298,360,369,341,310,260,369,371,279,357,349,296,307,349,284,288,348,305,338	ef:Z:ST-E00283:97:HY3WCCCXX:8:2214:24860:3173 ST-E00283:97:HY3WCCCXX:5:2205:31284:70733 ST-E00283:97:HY3WCCCXX:1:1213:9851:2979 ST-E00283:97:HY3WCCCXX:8:2207:1397:4684 ST-E00283:97:HY3WCCCXX:5:1120:6989:24075 ST-E00283:97:HY3WCCCXX:2:1222:27001:41515 ST-E00283:97:HY3WCCCXX:3:1210:14600:28154 ST-E00283:97:HY3WCCCXX:8:1218:8948:1836 ST-E00283:97:HY3WCCCXX:6:2122:27712:4772 ST-E00283:97:HY3WCCCXX:7:2213:8166:62786 ST-E00283:97:HY3WCCCXX:2:2112:14793:8798 ST-E00283:97:HY3WCCCXX:7:2222:6309:10486 ST-E00283:97:HY3WCCCXX:2:1214:9414:29015 ST-E00283:97:HY3WCCCXX:2:2218:12997:17781 ST-E00283:97:HY3WCCCXX:4:1212:21734:26800 ST-E00283:97:HY3WCCCXX:2:2219:11769:68394 ST-E00283:97:HY3WCCCXX:2:1210:26595:8851 ST-E00283:97:HY3WCCCXX:5:2107:1763:22440 ST-E00283:97:HY3WCCCXX:2:1209:19390:52906 ST-E00283:97:HY3WCCCXX:6:2224:12672:28294 ST-E00283:97:HY3WCCCXX:4:2104:13738:69519 ST-E00283:97:HY3WCCCXX:8:1111:24637:59499 ST-E00283:97:HY3WCCCXX:7:1102:14103:34535 ST-E00283:97:HY3WCCCXX:7:2221:5223:28892 ST-E00283:97:HY3WCCCXX:6:2207:30523:70820 ST-E00283:97:HY3WCCCXX:6:1108:9922:54436 ST-E00283:97:HY3WCCCXX:2:1207:31395:47334 ST-E00283:97:HY3WCCCXX:5:2102:25865:41620 ST-E00283:97:HY3WCCCXX:2:1214:20141:9413 ST-E00283:97:HY3WCCCXX:1:1209:5233:45576 ST-E00283:97:HY3WCCCXX:6:1217:18710:11998 ST-E00283:97:HY3WCCCXX:3:1106:17817:2118 ST-E00283:97:HY3WCCCXX:4:1103:4696:6460 ST-E00283:97:HY3WCCCXX:2:2211:19461:60765 ST-E00283:97:HY3WCCCXX:7:1217:12773:62470 ST-E00283:97:HY3WCCCXX:4:1112:9272:72192 ST-E00283:97:HY3WCCCXX:2:2209:8877:29595 ST-E00283:97:HY3WCCCXX:7:1116:25002:64667 ST-E00283:97:HY3WCCCXX:2:1215:1580:61679 ST-E00283:97:HY3WCCCXX:1:1216:22282:33516 ST-E00283:97:HY3WCCCXX:3:2216:11525:43431 ST-E00283:97:HY3WCCCXX:4:1104:31182:4667 ST-E00283:97:HY3WCCCXX:4:1222:25256:27134 ST-E00283:97:HY3WCCCXX:3:2108:13859:65195 ST-E00283:97:HY3WCCCXX:2:2221:17229:19487 ST-E00283:97:HY3WCCCXX:8:1212:31040:2417 ST-E00283:97:HY3WCCCXX:7:2204:13027:47193 ST-E00283:97:HY3WCCCXX:6:1205:24840:18221 ST-E00283:97:HY3WCCCXX:3:2118:27722:66390 ST-E00283:97:HY3WCCCXX:1:2201:18893:28664 ST-E00283:97:HY3WCCCXX:1:2120:12266:72561 ST-E00283:97:HY3WCCCXX:1:2224:14773:47298 ST-E00283:97:HY3WCCCXX:7:1110:26656:58462 ST-E00283:97:HY3WCCCXX:8:1123:7141:72227 ST-E00283:97:HY3WCCCXX:6:1223:18375:28540 ST-E00283:97:HY3WCCCXX:2:1118:12449:59833 ST-E00283:97:HY3WCCCXX:7:1209:2098:68236 ST-E00283:97:HY3WCCCXX:4:2205:29518:66759 ST-E00283:97:HY3WCCCXX:2:2117:26616:19118 ST-E00283:97:HY3WCCCXX:3:2102:11008:1784 ST-E00283:97:HY3WCCCXX:4:2214:11150:73299 ST-E00283:97:HY3WCCCXX:4:1116:27336:29261 ST-E00283:97:HY3WCCCXX:8:2215:24799:14740 ST-E00283:97:HY3WCCCXX:8:2216:10318:56141 ST-E00283:97:HY3WCCCXX:7:2220:17787:55051 ST-E00283:97:HY3WCCCXX:7:1220:14976:25569 ST-E00283:97:HY3WCCCXX:8:1119:21542:29806 ST-E00283:97:HY3WCCCXX:7:1223:31690:20032 ST-E00283:97:HY3WCCCXX:4:2201:22343:59429 ST-E00283:97:HY3WCCCXX:8:1222:20983:62487 ST-E00283:97:HY3WCCCXX:7:2218:22019:12666 ST-E00283:97:HY3WCCCXX:3:1113:6695:53100 ST-E00283:97:HY3WCCCXX:2:2203:12063:64122 ST-E00283:97:HY3WCCCXX:1:2218:31700:18573 ST-E00283:97:HY3WCCCXX:6:1219:29122:60905 ST-E00283:97:HY3WCCCXX:6:1113:27925:72438 ST-E00283:97:HY3WCCCXX:3:2217:8298:54647 ST-E00283:97:HY3WCCCXX:8:2213:15727:55245 ST-E00283:97:HY3WCCCXX:3:2122:9770:66759 ST-E00283:97:HY3WCCCXX:2:1207:9283:55825 ST-E00283:97:HY3WCCCXX:2:1224:28848:22071 ST-E00283:97:HY3WCCCXX:3:1102:1610:63454 ST-E00283:97:HY3WCCCXX:5:1213:27245:48125 ST-E00283:97:HY3WCCCXX:4:1223:15199:62909 ST-E00283:97:HY3WCCCXX:8:2220:2402:35783 ST-E00283:97:HY3WCCCXX:5:1107:10927:39757 ST-E00283:97:HY3WCCCXX:7:2104:20760:36012 ST-E00283:97:HY3WCCCXX:2:1117:26159:64316 ST-E00283:97:HY3WCCCXX:4:2120:22221:1590 ST-E00283:97:HY3WCCCXX:6:2215:27945:44028 ST-E00283:97:HY3WCCCXX:5:2112:24119:71612 ST-E00283:97:HY3WCCCXX:6:1106:29051:12191 ST-E00283:97:HY3WCCCXX:7:1215:27965:22651 ST-E00283:97:HY3WCCCXX:3:1202:15615:15320 ST-E00283:97:HY3WCCCXX:6:1110:12855:46261 ST-E00283:97:HY3WCCCXX:3:2108:29335:4175 ST-E00283:97:HY3WCCCXX:3:2108:29335:4175 ST-E00283:97:HY3WCCCXX:5:2206:21481:48371 ST-E00283:97:HY3WCCCXX:4:1106:1235:31230 ST-E00283:97:HY3WCCCXX:5:2111:16975:9624 ST-E00283:97:HY3WCCCXX:1:2104:2970:22985 ST-E00283:97:HY3WCCCXX:1:2116:6756:44275 ST-E00283:97:HY3WCCCXX:7:2216:29183:37207 ST-E00283:97:HY3WCCCXX:2:2112:8856:36557 ST-E00283:97:HY3WCCCXX:1:2201:5416:56405 ST-E00283:97:HY3WCCCXX:3:2207:22546:59112 ST-E00283:97:HY3WCCCXX:3:2113:8998:34904 ST-E00283:97:HY3WCCCXX:8:1219:23236:71524 ST-E00283:97:HY3WCCCXX:1:1104:11515:5792 ST-E00283:97:HY3WCCCXX:8:1214:23023:37436 ST-E00283:97:HY3WCCCXX:5:2210:5071:60589 ST-E00283:97:HY3WCCCXX:3:1121:9475:37313 ST-E00283:97:HY3WCCCXX:4:1123:15737:61099 ST-E00283:97:HY3WCCCXX:2:2107:19624:31265 ST-E00283:97:HY3WCCCXX:4:1109:30168:34061 ST-E00283:97:HY3WCCCXX:4:2114:11038:17307 ST-E00283:97:HY3WCCCXX:2:1116:19715:67779 ST-E00283:97:HY3WCCCXX:8:1214:11383:42517	dl:i:151	eq:B:f,23.2622,23.2622,23.2622,23.2622,23.2622,23.2622,23.2622,23.2622,23.0597,23.2622,23.0597,23.0597,23.2622,23.2622,23.2622,23.2622,23.0597,16.1564,23.2622,23.0597,23.2622,23.2622,23.2622,23.2622,23.0597,23.2622,23.2622,15.6654,15.817,23.2622,23.2622,23.0597,23.2622,23.2622,23.0597,23.2622,23.2622,23.0597,23.0597,23.2622,23.2622,23.0597,23.2622,23.0597,15.6654,23.0597,23.2622,23.2622,23.2622,23.0597,23.0597,23.0597,23.0597,23.2622,23.0597,23.0597,23.0597,23.2622,23.2622,23.0597,23.2622,23.0597,23.2622,15.9764,23.2622,23.0597,23.0597,23.0597,23.0597,23.2622,23.2622,23.2622,23.0597,23.0597,23.2622,23.2622,23.2622,23.0597,23.2622,23.2622,23.0597,23.2622,23.2622,23.2622,23.0597,23.2622,23.0597,23.2622,23.2622,23.2622,23.2622,23.2622,23.2622,23.2622,23.2622,17.4742,23.2622,23.2622,23.2622,23.0597,23.2622,23.0597,23.2622,23.2622,23.2622,23.0597,23.2622,23.2622,23.2622,23.0597,23.2622,23.2622,23.2622,23.2622,23.0597,23.0597,23.2622,23.0597	os:B:i,192,105,213,123,177,186,168,161,206,150,183,206,227,186,224,202,193,235,137,192,214,147,191,207,152,117,132,237,236,134,168,171,227,199,235,136,233,196,218,113,202,213,190,206,237,228,169,172,113,222,129,155,227,112,174,123,135,151,209,228,143,159,215,235,144,221,160,145,109,120,173,154,161,206,165,216,202,221,146,196,230,155,219,220,162,221,166,174,114,169,173,114,227,151,230,0,221,223,115,149,211,220,204,161,111,220,222,130,208,215,147,158,200,135,139,199,156,189	at:B:i,0,0	et:B:c,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0	ez:Z:-XuMsgbcmPuQZcp3japBqZn2rj_to7ZO -tfsOZrGGro76RxsImrf2CJp8jr9_ScC -wzH3_-x2yOCDTO8Rm6RmI0iy78qd5ZP 13j3zQWtHi3TuYCQhRX9zfjiREkeWMi5 1J208jgXgi1VOrKnGB_q6Iu9Kp3pQI37 1o_ulwO9eu3ubti4eYv6M2xQp38DiMgV 3TgfJaV2zEvYLSmrwyQWoK_YyXR45O3g 3qhnWMAsGm-vmR7EaiDMJkjve_5tkoRp 5VNBe43pLJzwXD8s4y3oOCplZ5dZIh0A 5ZXLUEgFMeWoxpGlIITHo030UhwLunkS 5nOhxvyBlrmRkwm83tl_7BC7oNRCaFzJ 6hCP7rz1Y8VdTYHkZ8wxGIzXX3aMchJK 6pk-tgXu7dBGBPOARNaw_e2BTas-b26D 7KtBFSU5eBc2pwioUZmnjiFMvQD_2JBh 9BPUSnmVEe0-kkhXL1nofm2B3J1ZmhRW 9bSztljh9SSUl_diVVRJoe0tbcFwzjkK A1waEMSNc6U9D4o8T79lQjfhU4HCeKSM ARKWfZz776SdiJ-C7gXi21Gz7Lb1I5gc AXyY_1u0EDOL1j9ksteEA9dk6A6GJ9ry AmGbrS-dfHFcXj55stBUyk03pL8YYOaK BMRSPDWV75PFJcdCGZZhZrkD0_TXZf2r Bvr5Hdk1dOnAtxe8Foq_zxmV04DAWmRJ CZ982hQowM2kz2kr96J2Gk-i6W4GXO_y Cl2-GFjECYxi5XvZTOPfi2Cku2eIiym2 CyCZcIxLqnKzy9nBtDJOHEJBX5LG_cqJ Du1L5elw2wJmfdLkOqT0LrUlvt-carxo EMb9NAE_U6HRa48hT3ftwtanW1eGZ6IR EiNd51Bt_CLdX7fbNYDh_3RV7LPsHXVA GFiegvUWcMqg8Nv2fFbltHUP9pFVGl4p GIcSbVn4XXzwEakrjqLztDBTsHh3qf0h GQUIcwjeWpsBMkM-_q82T-mRd7x2I-H4 HL1w6ku7PYtqvsOVZKyU8PP7WopaezD3 INw_6DuKR4A9GXMbylJSJQoAFC0BFjVx IPNLXKWUmT0lzA5ZHVMXmcWUB8j8SLDJ I_B6Kh8KBTOrgyXrNKGPqSIrkV5Bdy5g JQoRT8qhIEx5oUqpsOYt2XpZRG_rPkGu JkTNPSyvdX5ULAXOoWChsRkP-W_5JqB6 KBS-76M8FAxO9hJS3SCZI_Z6-BHKs4Ot KV7falKm1u8jDwd37hzNGzU3rE78CEqr MfvR-mvmTC2Hj3j9yUF7RH_G71CSHP7T MiKNE3jXdCrjjTocvDNPPbdgjUC-S7pg MyjYliy_sI457cx8DdC2V8WqkTMIRCSD NszBKESkBda4-mbmmheNheTWPNJ3hngQ Ok0ja2-ByOphbbHjjQNUYrFq-Z3MeWDq OoGH3h1-MMOsMQOxWKSE4Q6PGKpjC1uc PAfNFSN083TiIK3stol5II9xX6UWYsHo QVT_SNZpYEmg6GFHaH8OILGWZdopB0fr Qgl7GZj6PpqDJtZWWTJNXUUS2-oxEOja Qu-mHUa7ar37LPU0JUt3l9Xq6ebBqXvd RUYx7UiFAjlGJljrbnWJ2lZ4YiCIIeDq R_4II1Qnj0iY-mpr4AADkP9T-BBPaD2O RxSABSPe0_9-_bti2kuWETvxkI_5iVbS SVzqmxvjKkq4fCOyOOaFq9uIJU6V9CSy SaFSOqibhBXoCsipUDWOuGIlYUAIWSNA T7pwbMea7_wYyM4o3a3o-7nX7ViKB0lF TUwmW-61xthAmFfKhcbI_JA_CorIYh0U Tqb3pJawHFyycmsf-0qAVc2B3RCiX-sz TunOVna29SJCRUZlhYtYDRFnsISBJmBf Tw2KLEAK0P1kI7wDUf6irt5qQtA-6SUA TwWKUCBfNokcUh43yUyqkJRHnVtZrDdD W-YRZ23sUrUHfwoQGlFitB6Z0GEa-0Ro W8b4pbL6eH-weZScf5gftgyR1XXuoBcc X8mtwn6exm5gKBcuLrSQDWNnzlGNoCeO XpUKI68hfUP-DeaW9la7qRmk-19fSir1 Y_7N5V7gUVc_QNsb_ur9AoU8hKiQEMnd _Dz_XSQ5rQzwSkue956bou0r2WOot5od _NOlHDS4b9d-l8v5PRK-7z2BAeaqNv6K _Nvf2WW1id351gKyaDCavhfyVH1fPCMr aFPeIAyV77tn0FULH5SW5Ou6KIXZQicc aGD1Op8A4Icf9bSckSR3Z8WA2J1b7-Pg aT2BFhVE2lgs6D_GvshSS5P3SlN_zNbu b868rX_24-RMf44YqbR3WOwbFLS_Zlms b8hEgPcO3pyVmBNckg6m3JXIgJ-nHG8Y bEe5_G7yWsU6ETroRNUULqM-CtO0hbmQ cUY2c68KJiikl8u-abCyjh5ZO37Daa80 chQvn8ePX-LssULNsHb1yO5IM8KwNxjf dktum6PnlqWaJxFAq0CkXl0Pwyu4aJ7T dvY5XdRsc3yWqsQS09JsP2GS5rAj6FVh fItTw2YKMJy9Pgaj4DCmT4UWdi_UQA7J fPSbnjzuIgwfS1bbWal2mpoupCnGIpXS fnhGzFXanGEuW73Zc-5gtw10FSfQ90NK ghNjt6dy4zOKcvK-VI8esq3xI5CM5XHv hJFWQ0KwKqUhUVvDzdIr4fPPku9_Ua1m hUj3RobFPZop_DZnwPuzYiayQKr__7YC jaJUcjtIKkljMBEHwD20RyWyaldMpyuQ kdKZzcrvrGsdkKh4-aPWpBwgHgIo3WaC lPGktKZ-VeuOf84_2SCh1L862vTYK9gn m2Zxp9ZkgSozXt_YXJuydOzwboWslm1j mQ3JDvzsKW3cCD0rAgFEgTZ526Sa3FH- mTKxZSA44L25Yh4Ket5zjqb-YAjMIEh8 mvTDsmS-BovCth4iv6w6rlF3KydAx3no n3Cc6kgBlk_VsGQja4LtOdeauWU0BIyL nLXO4pjaDCD4eEgdk0xII0no7mNy3QA7 nXJHryxqx3dSUxFVSTFrkSzc5AcZ-rmf o6Ln5wlv_7_wbstLwvgiYLxbL7VIae1C p166uAvMOKCcvxzhElXHsPC6TAF3fyhi p166uAvMOKCcvxzhElXHsPC6TAbSk2Mr pL_0C9X255nMTUIEE8ryPed3Miv69Zlg pNR5dQ9ZyiDzYccEqsqqZYpJAJZh5B7j pNfsTJ7eLJfjayz7QQTvNuDhgAoh-FZo pUGSUadU_Zri9rw4-yRkn8EbYbaL4Ria phkWq0CjKySe6_PeP1tDiVsOSch16ZqX rawkoFfPxcgff8IiMnzV0YPmLQXS0BD5 s75klgIfhIpyngAzwqVN4jOjgcG2IDa9 sMaWItdTBypd3Bin6j8NCoIMzaVMts_B tAyZ8RKgGPQJZ5qqvP420_cYqzOlWumZ uTdkracCWrzZgYWbJ-Q38JOkdU5rvh9S v8m2Qjgr-svbY3wDAf_llgDuvmYNFWbL vVQno251d4wSiJgrzFR4Mu0rWhUlvSRk vtk9vYc70a8dNBiwTwpoOceRI18ePbwF vvHiP07m-lDf3WeG-ihuvUcYRzIclx9m wzEORUXnCKua_Vzr9PMkYZmj2Z1lbqHg x6wFEEY4DLG8iDHTu4dS0BQdGjok2j9T xNMvrGIGTXdqkIgtC2XYlp75ynFsTxRp y8PXKXSFL2Fis-Z-mkXNSglBNm9l02Es yNvU5gOR0hypSwJMnVpW5FWBYhTvx4NW yvPbXnOXIjW-CER1hwsQpjTljfYXrhcy zccgFBGZckgOdRqI4_Dj8FFa7EALLSo_";
        final SAMLineParser parser = new SAMLineParser(new DefaultSAMRecordFactory(), ValidationStringency.DEFAULT_STRINGENCY, getHeader(), null, null);
        SAMRecord r = parser.parseLine(primary);
        AssemblyAttributes aa = new AssemblyAttributes(r);
        List<SplitReadEvidence> e = SplitReadEvidence.create(SES(), r);
        assertEquals(1, e.size());
        double qual = e.get(0).getBreakpointQual();
        assertTrue(qual < 100);
    }
}
