package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;

import au.edu.wehi.idsv.picard.InMemoryReferenceSequenceFile;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class AssemblyEvidenceSourceTest extends IntermediateFilesTest {
	private File assemblyFile;
	@Before
	public void setup() throws IOException {
		super.setup();
		assemblyFile = new File(super.testFolder.getRoot(), "breakend.bam");
	}
	@Test
	public void should_write_breakend_bam() throws IOException {
		createInput(RP(0, 1, 2, 1));
		SAMEvidenceSource ses = new SAMEvidenceSource(getCommandlineContext(), input, null, 0);
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(getCommandlineContext(), ImmutableList.of(ses), assemblyFile);
		aes.assembleBreakends(null);
		assertTrue(assemblyFile.exists());
	}
	@Test
	public void breakend_bam_should_be_coordinate_sorted() throws IOException {
		createInput(RP(0, 1, 2, 1));
		SAMEvidenceSource ses = new SAMEvidenceSource(getCommandlineContext(), input, null, 0);
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(getCommandlineContext(), ImmutableList.of(ses), assemblyFile);
		aes.assembleBreakends(null);
		try (SamReader r = SamReaderFactory.make().open(assemblyFile)) {
			assertEquals(SortOrder.coordinate, r.getFileHeader().getSortOrder());
		}
	}
	@Test
	public void debruijn_should_generate_bam() throws IOException {
		createInput(
				withSequence("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGA", Read(0, 1, "41M58S")),
				withSequence("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGAT", Read(0, 1, "41M59S"))
				);
		ProcessingContext pc = getCommandlineContext();
		pc.getConfig().getAssembly().minReads = 1;
		SAMEvidenceSource ses = new SAMEvidenceSource(pc, input, null, 0);
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, ImmutableList.of(ses), assemblyFile);
		aes.assembleBreakends(null);
		
		assertEquals(1, getRecords(assemblyFile).size());
		assertEquals("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGAT", S(getRecords(assemblyFile).get(0).getReadBases()));
	}
	@Test
	public void iterator_should_return_in_chr_order() throws IOException {
		createInput(
				withSequence("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGA", Read(0, 93, "41M58S")),
				withSequence("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGAT", Read(0, 93, "41M59S")),
				withSequence("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGA", Read(1, 95, "41M58S")),
				withSequence("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGAT", Read(1, 95, "41M59S")),
				withSequence("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGA", Read(2, 97, "41M58S")),
				withSequence("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGAT", Read(2, 97, "41M59S")),
				withSequence("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGA", Read(3, 99, "41M58S")),
				withSequence("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGAT", Read(3, 99, "41M59S"))
				);
		ProcessingContext pc = getCommandlineContext();
		pc.getConfig().getAssembly().minReads = 1;
		//pc.getRealignmentParameters().requireRealignment = false;
		SAMEvidenceSource ses = new SAMEvidenceSource(pc, input, null, 0);
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, ImmutableList.of(ses), assemblyFile);
		aes.assembleBreakends(null);
		List<DirectedEvidence> list = Lists.newArrayList(aes.iterator());
		assertEquals(4, list.size());
		for (int i = 0; i <= 3; i++) {
			assertEquals(i, list.get(i).getBreakendSummary().referenceIndex);
		}
	}
	@Test
	public void should_not_write_filtered_assemblies() throws IOException {
		createInput(
				OEA(0, 1, "100M", true)
				);
		ProcessingContext pc = getCommandlineContext();
		pc.getAssemblyParameters().writeFiltered = false;
		SAMEvidenceSource ses = new SAMEvidenceSource(pc, input, null, 0);
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, ImmutableList.of(ses), assemblyFile);
		aes.assembleBreakends(null);
		List<DirectedEvidence> contigs = Lists.newArrayList(aes.iterator());
		assertEquals(0, contigs.size());
	}
	@Test
	public void should_filter_fully_reference_assemblies() {
		SAMRecord r = AssemblyFactory.createAnchoredBreakend(
				getContext(), AES(), new SequentialIdGenerator("asm"), BWD, null,
				0, 1, 2, B("AA"), B("AA"));
		assertTrue(AES().shouldFilterAssembly(r));
	}
	@Test
	public void should_filter_single_read_assemblies() {
		MockDirectedEvidence ev = new MockDirectedEvidence(new BreakendSummary(0, FWD, 1, 1, 2));
		SAMRecord r = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), FWD, Lists.newArrayList(ev), 0, 1, 1, B("AA"), B("AA"));
		assertTrue(AES().shouldFilterAssembly(r));
	}
	@Test
	public void should_filter_mate_anchored_assembly_shorter_than_read_length() {
		ArrayList<DirectedEvidence> support = Lists.<DirectedEvidence>newArrayList(NRRP(OEA(0, 1, "4M", true)));
		SAMRecord e = AssemblyFactory.createUnanchoredBreakend(
				getContext(), AES(), new SequentialIdGenerator("asm"), new BreakendSummary(0, FWD, 5, 5, 300), support,
				B("AAA"), B("AAA"), new int[] { 2, 0 });
		assertTrue(AES().shouldFilterAssembly(e));
	}
	@Test
	public void read_length_filter_should_not_apply_to_anchored_breakend_assembly() {
		ArrayList<DirectedEvidence> support = Lists.<DirectedEvidence>newArrayList(SCE(FWD, Read(0, 1, "1M2S")),
				NRRP(OEA(0, 1, "3M", true)),
				NRRP(OEA(0, 1, "4M", true)));
		SAMRecord e = AssemblyFactory.createAnchoredBreakend(
				getContext(), AES(), new SequentialIdGenerator("asm"), FWD, support,
				0, 1, 1, B("AAA"), B("AAA"));
		assertFalse(AES().shouldFilterAssembly(e));
	}
	@Test
	public void soft_clip_size_filter_should_not_apply_to_unanchored_assembly() {
		ArrayList<DirectedEvidence> support = Lists.<DirectedEvidence>newArrayList(
				NRRP(OEA(0, 1, "3M", false)),
				NRRP(OEA(0, 1, "4M", false)));
		SAMRecord e = AssemblyFactory.createUnanchoredBreakend(
				getContext(), AES(), new SequentialIdGenerator("asm"), new BreakendSummary(0, FWD, 5, 5, 300), support,
				B("AAAAAA"), B("AAAAAA"), new int[] { 2, 0 });
		assertFalse(AES().shouldFilterAssembly(e));
	}
	@Test
	public void should_filter_if_no_breakpoint_assembly() {
		ArrayList<DirectedEvidence> support = Lists.<DirectedEvidence>newArrayList(
				SCE(BWD, Read(0, 1, "1S1M")),
				NRRP(OEA(0, 1, "3M", false)),
				NRRP(OEA(0, 1, "4M", false)));
		SAMRecord e = AssemblyFactory.createAnchoredBreakend(
				getContext(), AES(), new SequentialIdGenerator("asm"), BWD, support,
				0, 1, 2, B("AA"), B("AA"));
		// reference assembly
		assertTrue(AES().shouldFilterAssembly(e));
	}
	@Test
	public void should_not_apply_breakend_filter_to_unanchored_assembly() {
		AssemblyEvidenceSource aes = AES();
		aes.getContext().getConfig().getAssembly().minReads = 0;
		ArrayList<DirectedEvidence> support = Lists.<DirectedEvidence>newArrayList(NRRP(SES(100, 100), DP(0, 1, "1M", true, 0, 5, "1M", false)));
		SAMRecord e = AssemblyFactory.createUnanchoredBreakend(
				getContext(), aes, new SequentialIdGenerator("asm"), new BreakendSummary(0, FWD, 1, 1, 300), support,
				B("AA"), B("AA"), new int[] { 2, 0});
		
		assertFalse(aes.shouldFilterAssembly(e));
	}
	@Test
	public void should_filter_too_few_reads() {
		ProcessingContext pc = getContext();
		AssemblyEvidenceSource aes = AES(pc);
		pc.getAssemblyParameters().minReads = 3;
		SAMRecord e = AssemblyFactory.createAnchoredBreakend(pc, aes, new SequentialIdGenerator("asm"), BWD, null, 0, 1, 5, B("AACGTG"), B("AACGTG"));
		assertTrue(aes.shouldFilterAssembly(e));
		
		ArrayList<DirectedEvidence> support = Lists.<DirectedEvidence>newArrayList(
				SCE(BreakendDirection.Forward, withQual(new byte[] { 5,5,5,5,5,5 }, withSequence("AACGTG", Read(0, 1, "1M5S")))));
		e = AssemblyFactory.createAnchoredBreakend(
				getContext(), aes, new SequentialIdGenerator("asm"), BWD, support,
				0, 1, 5, B("AACGTG"), B("AACGTG"));
		assertTrue(aes.shouldFilterAssembly(e));
		
		support = Lists.<DirectedEvidence>newArrayList(
				SCE(BreakendDirection.Forward, withQual(new byte[] { 5,5,5,5,5,5 }, withSequence("AACGTG", Read(0, 1, "1M5S")))),
				NRRP(OEA(0, 1, "4M", false)));
		e = AssemblyFactory.createAnchoredBreakend(
				pc, aes, new SequentialIdGenerator("asm"), BWD, support,
				0, 1, 5, B("AACGTG"), B("AACGTG"));
		assertTrue(aes.shouldFilterAssembly(e));
		
		support = Lists.<DirectedEvidence>newArrayList(
				SCE(BreakendDirection.Forward, withQual(new byte[] { 5,5,5,5,5,5 }, withSequence("AACGTG", Read(0, 1, "1M5S")))),
				NRRP(OEA(0, 1, "3M", false)),
				NRRP(OEA(0, 1, "5M", false)));
		e = AssemblyFactory.createAnchoredBreakend(
				pc, aes, new SequentialIdGenerator("asm"), BWD, support,
				0, 1, 5, B("AACGTG"), B("AACGTG"));
		assertFalse(AES().shouldFilterAssembly(e));
	}
	@Test
	public void should_filter_reference_breakend() {
		List<DirectedEvidence> evidence = Lists.newArrayList();
		evidence.add(SCE(FWD, Read(0, 5, "5M5S")));
		//evidence.add(SCE(FWD, Read(0, 5, "5M6S")));
		//evidence.add(SCE(FWD, Read(0, 5, "5M7S")));
		SAMRecord e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), FWD, evidence,
				0, 10, 5, B("AAAAAAAAAA"), B("AAAAAAAAAA"));
		assertTrue(AES().shouldFilterAssembly(e));
	}
	@Test
	public void should_filter_breakend_matching_reference_due_to_homology() {
		List<DirectedEvidence> evidence = Lists.newArrayList();
		evidence.add(SCE(FWD, Read(0, 5, "5M5S")));
		evidence.add(SCE(FWD, Read(0, 5, "5M6S")));
		evidence.add(SCE(FWD, Read(0, 5, "5M7S")));
		SAMRecord e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), FWD, evidence,
				0, 10, 5, B("AAAAAAAAAA"), B("AAAAAAAAAA"));
		SAMRecordUtil.unclipExactReferenceMatches(SMALL_FA, e);
		assertTrue(AES().shouldFilterAssembly(e));
	}
	@Test
	public void should_adjust_assembly_contig_alignment_to_align_with_reference_homology() throws IOException {
		createInput(
				Read(0, 5, "5M5S"),
				Read(0, 5, "5M6S"),
				Read(0, 5, "5M7S")
				);
		SAMEvidenceSource ses = new SAMEvidenceSource(getCommandlineContext(), input, null, 0);
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(getCommandlineContext(), ImmutableList.of(ses), assemblyFile);
		aes.assembleBreakends(null);
		List<SAMRecord> assemblies = getRecords(assemblyFile);
		assertEquals(0, assemblies.size());
	}
	@Test
	@Ignore("Enhancement") // TODO: should we even do this? How much of this is covered by SingleReadEvidence.isReference()?
	public void should_realign_very_poor_matches_to_the_reference() {
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
		ProcessingContext pc = new ProcessingContext(new FileSystemContext(testFolder.getRoot(), 500000), reference, ref, Lists.newArrayList(), getConfig(testFolder.getRoot()));
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, ImmutableList.of(SES(pc)), assemblyFile);
		SAMRecord ra = aes.transformAssembly(r);
		assertEquals(2, ra.getAlignmentStart());
		assertEquals("131M", ra.getCigarString());
	}
	@Test
	public void chunck_spanning_assemblies_should_not_be_repeated() throws IOException {
		List<SAMRecord> in = new ArrayList<>();
		for (int i = 50; i < 150; i++) {
			in.add(withSequence("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGA", Read(0, i, "41M58S"))[0]);
		}
		createInput(in);
		ProcessingContext pc = getCommandlineContext();
		pc.getConfig().getAssembly().minReads = 1;
		pc.getConfig().chunkSize = 100;
		SAMEvidenceSource ses = new SAMEvidenceSource(pc, input, null, 0);
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, ImmutableList.of(ses), assemblyFile);
		aes.assembleBreakends(null);
		List<DirectedEvidence> list = Lists.newArrayList(aes.iterator());
		assertEquals(100, list.size());
	}
	@Test
	public void parallel_assembly_should_not_affect_assembly_results() throws IOException {
		List<SAMRecord> in = new ArrayList<>();
		for (int i = 50; i < 150; i++) {
			in.add(withSequence("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGA", Read(0, i, "41M58S"))[0]);
		}
		createInput(in);
		ProcessingContext pc = getCommandlineContext();
		pc.getConfig().getAssembly().minReads = 1;
		pc.getConfig().chunkSize = 100;
		SAMEvidenceSource ses = new SAMEvidenceSource(pc, input, null, 0);
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, ImmutableList.of(ses), assemblyFile);
		ExecutorService threadpool = Executors.newFixedThreadPool(8);
		aes.assembleBreakends(threadpool);
		threadpool.shutdown();
		List<DirectedEvidence> list = Lists.newArrayList(aes.iterator());
		assertEquals(100, list.size());
	}
	@Test
	public void should_deduplicate_assembly_with_multimapping_reads() throws IOException {
		createInput(
				withReadName("read1", withSequence("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGA", Read(0, 1, "41M58S"))),
				withReadName("read2", withSequence("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGAT", Read(0, 1, "41M59S"))),
				applying(f -> f.setNotPrimaryAlignmentFlag(true), withReadName("read1", withSequence("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGA", Read(0, 1000, "41M58S")))),
				applying(f -> f.setNotPrimaryAlignmentFlag(true), withReadName("read2", withSequence("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGAT", Read(0, 1000, "41M59S"))))
				);
		ProcessingContext pc = getCommandlineContext();
		pc.getConfig().getAssembly().minReads = 1;
		SAMEvidenceSource ses = new SAMEvidenceSource(pc, input, null, 0);
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, ImmutableList.of(ses), assemblyFile);
		aes.assembleBreakends(null);
		assertEquals(1, getRecords(assemblyFile).size());
	}
}
