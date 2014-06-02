package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqRecord;

import java.util.List;

import org.junit.Test;

public class GenerateDirectedBreakpointsTest extends CommandLineTest {
	public SAMRecord validSC(String seq, int referenceIndex, int position, String cigar) {
		SAMRecord r = withSequence(seq, Read(referenceIndex, position, cigar))[0];
		r.setMappingQuality(5);
		return r;
	}
			
	public SAMRecord ValidSC() {
		SAMRecord r = Read(0, 1, "25S50M25S");
		r.setReadName("SC1");
		byte[] qual = new byte[100];
		for (int i = 0; i < qual.length; i++) {
			qual[i] = 5;
		}
		r.setBaseQualities(qual);
		r.setMappingQuality(5);
		return r;
	}
	@Test
	public void should_write_fastq() {
		createInput(RP(0, 1, 2, 1));
		extractEvidence();
		generateDirectedBreakpoints("polyA");
		assertTrue(fastq.exists());
	}
	@Test
	public void should_write_long_sc_to_fastq() {
		createInput(ValidSC());
		extractEvidence();
		generateDirectedBreakpoints("polyA");
		assertTrue(fastq.exists());
		assertEquals(2, getFastqRecords().size());
		assertTrue(getFastqRecords().get(0).getReadHeader().contains("SC1"));
		assertTrue(getFastqRecords().get(1).getReadHeader().contains("SC1"));
		assertEquals("AAAAAAAAAAAAAAAAAAAAAAAAA", getFastqRecords().get(0).getReadString());
		// +33 encoding of mapping quality 5
		assertEquals("&&&&&&&&&&&&&&&&&&&&&&&&&", getFastqRecords().get(0).getBaseQualityString());
	}
	@Test
	public void min_mapq_should_filter_sc() {
		createInput(ValidSC());
		extractEvidence();
		generateDirectedBreakpoints("polyA", 6, null, null, null, null);
		assertTrue(fastq.exists());
		assertEquals(0, getFastqRecords().size());
	}
	@Test
	public void long_sc_length_should_filter_sc() {
		createInput(ValidSC());
		extractEvidence();
		generateDirectedBreakpoints("polyA", null, 26, null, null, null);
		assertTrue(fastq.exists());
		assertEquals(0, getFastqRecords().size());
	}
	@Test
	public void min_id_should_filter_sc() {
		SAMRecord r = ValidSC();
		byte[] readBases = r.getReadBases();
		for (int i = 25; i < 51; i++) {
			// change just over half the bases from A to T
			readBases[i] = 'T';
		}
		r.setReadBases(readBases);
		createInput(RP(0, 1, 2, 1));
		extractEvidence();
		generateDirectedBreakpoints("polyA", null, null, 50f, null, null);
		assertTrue(fastq.exists());
		assertEquals(0, getFastqRecords().size());
	}
	@Test
	public void min_qual_should_filter_sc() {
		createInput(ValidSC());
		extractEvidence();
		generateDirectedBreakpoints("polyA", null, null, null, 6f, null);
		assertTrue(fastq.exists());
		assertEquals(0, getFastqRecords().size());
	}
	@Test
	public void min_qual_should_only_consider_sc_qualities() {
		SAMRecord r = ValidSC();
		byte[] qual = r.getBaseQualities();
		for (int i = 25; i < 75; i++) {
			qual[i] = 40;
		}
		r.setBaseQualities(qual);
		createInput(r);
		extractEvidence();
		generateDirectedBreakpoints("polyA", null, null, null, 6f, null);
		assertTrue(fastq.exists());
		assertEquals(0, getFastqRecords().size());
	}
	@Test
	public void debruijn_should_generate_fastq() {
		createInput(
				withSequence("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGA", Read(0, 1, "1M98S")),
				withSequence("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGAT", Read(0, 1, "1M99S"))
				);
		extractEvidence();
		generateDirectedBreakpoints("polyA", null, null, null, null, null);
		assertTrue(fastq.exists());
		assertEquals(1, getFastqRecords().size());
		assertEquals("ATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGAT", getFastqRecords().get(0).getReadString());
	}
	@Test
	public void debruijn_should_generate_vcf() {
		SAMRecord r1 = Read(0, 1, "1M5S");
		r1.setReadBases(B("AACGT"));
		SAMRecord r2 = Read(0, 1, "1M6S");
		r2.setReadBases(B("AACGTG"));
		createInput(r1, r2);
		extractEvidence();
		generateDirectedBreakpoints("polyA", null, 1, null, null, 3);
		assertEquals(1, getVcf(".idsv.polyA.breakend.vcf").size());
	}
	@Test
	public void should_write_vcf() {
		createInput(RP(0, 1, 2, 1));
		extractEvidence();
		generateDirectedBreakpoints("polyA");
		workingFileShouldExist(".idsv.polyA.breakend.vcf");
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
		extractEvidence();
		generateDirectedBreakpoints("polyA");
		assertTrue(fastq.exists());
		List<FastqRecord> out = getFastqRecords();
		assertEquals(6, out.size());
		assertTrue(out.get(0).getReadHeader().startsWith("0#1#"));
		assertTrue(out.get(1).getReadHeader().startsWith("0#1#"));
		assertTrue(out.get(2).getReadHeader().startsWith("0#1#"));
		assertTrue(out.get(3).getReadHeader().startsWith("0#10#"));
		assertTrue(out.get(4).getReadHeader().startsWith("0#10#"));
		assertTrue(out.get(5).getReadHeader().startsWith("0#10#"));
	}
}
