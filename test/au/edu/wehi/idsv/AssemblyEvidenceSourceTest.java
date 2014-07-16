package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqRecord;

import java.io.File;
import java.util.List;

import org.junit.Test;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;

public class AssemblyEvidenceSourceTest extends IntermediateFilesTest {
	@Test
	public void should_write_fastq() {
		createInput(RP(0, 1, 2, 1));
		SAMEvidenceSource ses = new SAMEvidenceSource(getCommandlineContext(), input, false);
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(getCommandlineContext(), ImmutableList.of(ses), new File(super.testFolder.getRoot(), "out.vcf"));
		aes.ensureAssembled();
		getAssembly(aes);
		getFastqRecords(aes);
	}
	
	@Test
	public void debruijn_should_generate_fastq() {
		createInput(
				withSequence("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGA", Read(0, 1, "1M98S")),
				withSequence("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGAT", Read(0, 1, "1M99S"))
				);
		ProcessingContext pc = getCommandlineContext();
		pc.getAssemblyParameters().method = AssemblyMethod.DEBRUIJN_PER_POSITION;
		SAMEvidenceSource ses = new SAMEvidenceSource(pc, input, false);
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, ImmutableList.of(ses), new File(super.testFolder.getRoot(), "out.vcf"));
		aes.ensureAssembled();

		assertEquals(1, getFastqRecords(aes).size());
		assertEquals("ATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGAT", getFastqRecords(aes).get(0).getReadString());
	}
	@Test
	public void iterator_should_return_in_chr_order() {
		createInput(
				withSequence("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGA", Read(0, 1, "1M98S")),
				withSequence("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGAT", Read(0, 1, "1M99S")),
				withSequence("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGA", Read(1, 1, "1M98S")),
				withSequence("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGAT", Read(1, 1, "1M99S")),
				withSequence("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGA", Read(2, 1, "1M98S")),
				withSequence("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGAT", Read(2, 1, "1M99S")),
				withSequence("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGA", Read(3, 1, "1M98S")),
				withSequence("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGAT", Read(3, 1, "1M99S"))
				);
		ProcessingContext pc = getCommandlineContext();
		pc.getAssemblyParameters().method = AssemblyMethod.DEBRUIJN_PER_POSITION;
		SAMEvidenceSource ses = new SAMEvidenceSource(pc, input, false);
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, ImmutableList.of(ses), new File(super.testFolder.getRoot(), "out.vcf"));

		List<DirectedEvidence> list = Lists.newArrayList(aes.iterator());
		assertEquals(4,  list.size());
		for (int i = 0; i <= 3; i++) {
			assertEquals(i, list.get(i).getBreakendSummary().referenceIndex);
		}
	}
	@Test
	public void iterator_chr_should_return_only_chr_assemblies() {
		createInput(
				withSequence("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGA", Read(0, 1, "1M98S")),
				withSequence("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGAT", Read(0, 1, "1M99S")),
				withSequence("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGA", Read(1, 1, "1M98S")),
				withSequence("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGAT", Read(1, 1, "1M99S")),
				withSequence("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGA", Read(2, 1, "1M98S")),
				withSequence("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGAT", Read(2, 1, "1M99S")),
				withSequence("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGA", Read(3, 1, "1M98S")),
				withSequence("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGAT", Read(3, 1, "1M99S"))
				);
		ProcessingContext pc = getCommandlineContext();
		pc.getAssemblyParameters().method = AssemblyMethod.DEBRUIJN_PER_POSITION;
		SAMEvidenceSource ses = new SAMEvidenceSource(pc, input, false);
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, ImmutableList.of(ses), new File(super.testFolder.getRoot(), "out.vcf"));

		List<DirectedEvidence> list = Lists.newArrayList(aes.iterator("polyA"));
		assertEquals(1,  list.size());
	}
	@Test
	public void debruijn_should_generate_vcf() {
		SAMRecord r1 = Read(0, 1, "1M5S");
		r1.setReadBases(B("CAACGT"));
		SAMRecord r2 = Read(0, 1, "1M6S");
		r2.setReadBases(B("CAACGTG"));
		createInput(r1, r2);
		ProcessingContext pc = getCommandlineContext();
		pc.getAssemblyParameters().k = 3;
		pc.getAssemblyParameters().method = AssemblyMethod.DEBRUIJN_PER_POSITION;
		SAMEvidenceSource ses = new SAMEvidenceSource(pc, input, false);
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, ImmutableList.of(ses), new File(super.testFolder.getRoot(), "out.vcf"));
		aes.ensureAssembled();
		assertEquals(1, getAssembly(aes).size());
	}
	@Test
	public void should_write_vcf() {
		createInput(RP(0, 1, 2, 1));
		ProcessingContext pc = getCommandlineContext();
		pc.getAssemblyParameters().method = AssemblyMethod.DEBRUIJN_PER_POSITION;
		SAMEvidenceSource ses = new SAMEvidenceSource(pc, input, false);
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, ImmutableList.of(ses), new File(super.testFolder.getRoot(), "out.vcf"));
		aes.ensureAssembled();
		assertEquals(0, getAssembly(aes).size());
	}
	@Test
	public void should_make_calls_in_order() {
		createInput(
				// first assembly
				validSC("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGA", 0, 1, "1M98S"),
				validSC("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGAT", 0, 1, "1M99S"),
				// second assembly
				validSC("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGA", 0, 10, "1M98S"),
				validSC("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGAT", 0, 10, "1M99S")
				);
		ProcessingContext pc = getCommandlineContext();
		pc.getAssemblyParameters().method = AssemblyMethod.DEBRUIJN_PER_POSITION;
		SAMEvidenceSource ses = new SAMEvidenceSource(pc, input, false);
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, ImmutableList.of(ses), new File(super.testFolder.getRoot(), "out.vcf"));
		aes.ensureAssembled();
		assertEquals(2, getAssembly(aes).size());
		List<FastqRecord> out = getFastqRecords(aes);
		assertEquals(2, out.size());
		assertEquals(1, BreakpointFastqEncoding.getStartPosition(out.get(0).getReadHeader()));
		assertEquals(10, BreakpointFastqEncoding.getStartPosition(out.get(1).getReadHeader()));
	}
}
