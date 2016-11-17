package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.util.EnumSet;
import java.util.List;

import org.junit.Test;
import org.junit.experimental.categories.Category;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;

import au.edu.wehi.idsv.configuration.GridssConfiguration;
import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.picard.SynchronousReferenceLookupAdapter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.metrics.Header;
import htsjdk.samtools.metrics.StringHeader;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class AssemblyEvidenceSourceTest extends IntermediateFilesTest {
	@Test
	public void should_write_fastq() {
		createInput(RP(0, 1, 2, 1));
		SAMEvidenceSource ses = new SAMEvidenceSource(getCommandlineContext(), input, 0);
		ses.completeSteps(EnumSet.allOf(ProcessStep.class));
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(getCommandlineContext(), ImmutableList.of(ses), new File(super.testFolder.getRoot(), "out.vcf"));
		aes.ensureAssembled();
		getFastqRecords(aes);
	}
	@Test
	public void createUnanchoredBreakend_should_round_trip_unanchored_breakend_start_end_but_not_nominal_position() {
		assertEquals(new BreakendSummary(0, FWD, 2, 1, 3), AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new BreakendSummary(0, FWD, 1, 1, 3), null, B("GTAC"), new byte[] {1,2,3,4}, new int[] {0, 0}).getBreakendSummary());
		assertEquals(new BreakendSummary(0, BWD, 15, 10, 20), AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new BreakendSummary(0, BWD, 11, 10, 20), null, B("GTAC"), new byte[] {1,2,3,4}, new int[] {0, 0}).getBreakendSummary());
	}
	@Test
	public void debruijn_should_generate_fastq() {
		createInput(
				withSequence("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGA", Read(0, 1, "1M98S")),
				withSequence("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGAT", Read(0, 1, "1M99S"))
				);
		ProcessingContext pc = getCommandlineContext();
		SAMEvidenceSource ses = new SAMEvidenceSource(pc, input, 0);
		ses.completeSteps(ProcessStep.ALL_STEPS);
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, ImmutableList.of(ses), new File(super.testFolder.getRoot(), "out.vcf"));
		aes.ensureAssembled();

		assertEquals(1, getFastqRecords(aes).size());
		assertEquals("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGAT", getFastqRecords(aes).get(0).getReadString());
	}
	@Test
	public void iterator_should_return_in_chr_order() throws IOException {
		createInput(
				withSequence("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGA", Read(0, 1, "1M98S")),
				withSequence("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGAT", Read(0, 1, "1M99S")),
				withSequence("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGA", Read(1, 1, "1M98S")),
				withSequence("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGAT", Read(1, 1, "1M99S")),
				withSequence("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGA", Read(2, 1, "1M98S")),
				withSequence("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGAT", Read(2, 1, "1M99S")),
				withSequence("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGA", Read(3, 1, "1M98S")),
				withSequence("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGAT", Read(3, 1, "1M99S"))
				);
		ProcessingContext pc = getCommandlineContext();
		//pc.getRealignmentParameters().requireRealignment = false;
		SAMEvidenceSource ses = new SAMEvidenceSource(pc, input, 0);
		ses.completeSteps(ProcessStep.ALL_STEPS);
		writeUnmappedRealignments(pc, ses, 0);
		ses.completeSteps(ProcessStep.ALL_STEPS);
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, ImmutableList.of(ses), new File(super.testFolder.getRoot(), "out.vcf"));
		aes.ensureAssembled(); // assemble
		writeUnmappedRealignments(pc, aes, 0);
		aes.ensureAssembled(); // transform to annotated read pair

		List<SAMRecordAssemblyEvidence> list = Lists.newArrayList(aes.iterator(false, false));
		assertEquals(4,  list.size());
		for (int i = 0; i <= 3; i++) {
			assertEquals(i, list.get(i).getBreakendSummary().referenceIndex);
		}
	}
	@Test
	public void iterator_chr_should_return_only_chr_assemblies() throws IOException {
		createInput(
				withSequence("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGAC", Read(0, 1, "50M50S")),
				withSequence("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGAT", Read(0, 1, "49M51S")),
				withSequence("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGA", Read(1, 1, "50M49S")),
				withSequence("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGAT", Read(1, 1, "49M51S")),
				withSequence("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGA", Read(2, 1, "50M49S")),
				withSequence("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGAT", Read(2, 1, "49M51S")),
				withSequence("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGA", Read(3, 1, "50M49S")),
				withSequence("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGAT", Read(3, 1, "49M51S"))
				);
		ProcessingContext pc = getCommandlineContext();
		pc.getAssemblyParameters().subgraph.maxSubgraphFragmentWidth = -1;
		//pc.getRealignmentParameters().requireRealignment = false;
		SAMEvidenceSource ses = new SAMEvidenceSource(pc, input, 0);
		ses.completeSteps(ProcessStep.ALL_STEPS);
		writeUnmappedRealignments(pc, ses, 0);
		ses.completeSteps(ProcessStep.ALL_STEPS);
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, ImmutableList.of(ses), new File(super.testFolder.getRoot(), "out.vcf"));
		aes.ensureAssembled(); // assemble
		writeUnmappedRealignments(pc, aes, 0);
		aes.ensureAssembled(); // transform to annotated read pair

		List<SAMRecordAssemblyEvidence> list = Lists.newArrayList(aes.iterator(false, false, "polyA"));
		assertEquals(1,  list.size());
	}
	@Test
	public void debruijn_should_generate_vcf() {
		SAMRecord r1 = Read(0, 1, "1M5S");
		r1.setReadBases(B("AAACGT"));
		SAMRecord r2 = Read(0, 1, "1M6S");
		r2.setReadBases(B("AAACGTG"));
		createInput(r1, r2);
		ProcessingContext pc = getCommandlineContext();
		pc.getAssemblyParameters().k = 3;
		SAMEvidenceSource ses = new SAMEvidenceSource(pc, input, 0);
		ses.completeSteps(ProcessStep.ALL_STEPS);
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, ImmutableList.of(ses), new File(super.testFolder.getRoot(), "out.vcf"));
		aes.ensureAssembled();
		assertEquals(1, getAssemblyRaw(aes).size());
	}
	@Test
	public void should_write_raw_bam() {
		createInput(RP(0, 1, 2, 1));
		ProcessingContext pc = getCommandlineContext();
		SAMEvidenceSource ses = new SAMEvidenceSource(pc, input, 0);
		ses.completeSteps(EnumSet.allOf(ProcessStep.class));
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, ImmutableList.of(ses), new File(super.testFolder.getRoot(), "out.vcf"));
		aes.ensureAssembled();
		assertEquals(0, getAssemblyRaw(aes).size());
	}
	@Test
	public void should_make_calls_in_order() {
		createInput(
				// first assembly
				validSC("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGA", 0, 1, "1M98S"),
				validSC("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGAT", 0, 1, "1M99S"),
				// second assembly
				validSC("AATTAATCGCAAtAGCGGGAAGTATTCGACGCCCAGTCAGCTGGAGCACCATTACcCGATCAAtACATATCAGAtATGATTGACGTcTCACAAGCgGGA", 0, 10, "1M98S"),
				validSC("AATTAATCGCAAtAGCGGGAAGTATTCGACGCCCAGTCAGCTGGAGCACCATTACcCGATCAAtACATATCAGAtATGATTGACGTcTCACAAGCgGGAT", 0, 10, "1M99S")
				          //AGATCGGAAGAG
				);
		ProcessingContext pc = getCommandlineContext();
		SAMEvidenceSource ses = new SAMEvidenceSource(pc, input, 0);
		ses.completeSteps(ProcessStep.ALL_STEPS);
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, ImmutableList.of(ses), new File(super.testFolder.getRoot(), "out.vcf"));
		aes.ensureAssembled();
		assertEquals(2, getAssemblyRaw(aes).size());
		List<FastqRecord> out = getFastqRecords(aes);
		assertEquals(2, out.size());
		assertEquals(1, BreakpointFastqEncoding.getEncodedStartPosition(out.get(0).getReadHeader()));
		assertEquals(10, BreakpointFastqEncoding.getEncodedStartPosition(out.get(1).getReadHeader()));
	}
	@Test
	public void should_not_write_filtered_assemblies() {
		createInput(
				OEA(0, 1, "100M", true)
				);
		ProcessingContext pc = getCommandlineContext();
		pc.getAssemblyParameters().writeFiltered = false;
		SAMEvidenceSource ses = new SAMEvidenceSource(pc, input, 0);
		ses.completeSteps(ProcessStep.ALL_STEPS);
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, ImmutableList.of(ses), new File(super.testFolder.getRoot(), "out.vcf"));
		aes.ensureAssembled();
		assertEquals(0, getAssemblyRaw(aes).size());
	}
	@Test
	@Category(Hg19Tests.class)
	public void should_generate_indel_assemblies() throws IOException {
		createInput(new File("src/test/resources/inss.bam"));
		List<Header> headers = Lists.newArrayList();
		headers.add(new StringHeader("TestHeader"));
		File ref = Hg19Tests.findHg19Reference();
		IndexedFastaSequenceFile indexed = new IndexedFastaSequenceFile(ref);
		ReferenceLookup lookup = new SynchronousReferenceLookupAdapter(indexed);
		ProcessingContext pc = new ProcessingContext(
				new FileSystemContext(testFolder.getRoot(), 500000), ref, false, lookup,
				headers,
				new GridssConfiguration(getDefaultConfig(), testFolder.getRoot()));
		pc.getAssemblyParameters().anchorRealignment.perform = true;
		pc.registerCategory(0, "Normal");
		SAMEvidenceSource ses = new SAMEvidenceSource(pc, input, 0);
		ses.completeSteps(ProcessStep.ALL_STEPS);		
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, ImmutableList.of(ses), new File(super.testFolder.getRoot(), "out.vcf"));
		aes.ensureAssembled();
		assertEquals(2, getAssemblyRaw(aes).size());
		List<FastqRecord> out = getFastqRecords(aes);
		assertEquals(0, out.size());
		
	}
}
