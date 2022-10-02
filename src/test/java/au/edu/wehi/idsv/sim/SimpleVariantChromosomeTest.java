package au.edu.wehi.idsv.sim;

import au.edu.wehi.idsv.BreakpointSummary;
import au.edu.wehi.idsv.IdsvVariantContext;
import au.edu.wehi.idsv.IntermediateFilesTest;
import au.edu.wehi.idsv.bed.BedpeIterator;
import au.edu.wehi.idsv.bed.BedpeRecord;
import au.edu.wehi.idsv.vcf.SvType;
import com.google.common.collect.Lists;
import com.google.common.io.Files;
import htsjdk.samtools.SAMSequenceDictionary;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;


public class SimpleVariantChromosomeTest extends IntermediateFilesTest {
	public File outputBedpe;
	@Before
	public void setup() throws IOException {
		super.setup();
		outputBedpe = new File(testFolder.getRoot(), "SimpleVariantChromosomeTest.outputBedpe");
	}
	@After
	public void cleanup() throws Exception {
		super.cleanup();
		outputBedpe.delete();
	}

	@Test
	public void should_space_variants_with_margin_padding() throws IOException {
		SimpleVariantChromosome svc = new SimpleVariantChromosome(getContext(), "polyACGT", 1, 0, true);
		svc.assemble(input, output, outputBedpe, false, Lists.newArrayList(SvType.DEL, SvType.INS, SvType.INV, SvType.DUP), Lists.newArrayList(2), 1, false);
		// 123456  78901234  5
		// ACGTAC  GTACGTAC  G
		// ||  ||  ||  ||    |
		// AC--ACNNGTGTGTACACG
		//  *|| *|| *|| *||||
		//   v1  v2  v3   v4
		List<IdsvVariantContext> vcf = super.getVcf(output, null);
		assertEquals(2, vcf.get(0).getStart());
		assertEquals(6, vcf.get(1).getStart());
		assertEquals(8, vcf.get(2).getStart());
		assertEquals(12, vcf.get(3).getStart());
	}
	@Test
	public void should_output_breakpoint_positions() throws IOException {
		SimpleVariantChromosome svc = new SimpleVariantChromosome(getContext(), "polyACGT", 1, 0, true);
		SAMSequenceDictionary dict = svc.assemble(input, output, outputBedpe, false, Lists.newArrayList(SvType.DEL, SvType.INS, SvType.INV, SvType.DUP), Lists.newArrayList(2), 1, false);
		// 123456  78901234  5
		// ACGTAC  GTACGTAC  G
		// ||  ||  ||  ||    |
		// AC--ACNNGTGTGTACACG
		//  *|| *|| *|| *||||
		//   v1  v2  v3   v4
		// 12  345678901234567
		List<BedpeRecord> bps = Lists.newArrayList(new BedpeIterator(outputBedpe, dict));
		assertEquals(bps.get(0).bp, new BreakpointSummary(0, FWD, 2, 0, BWD, 3)); // DEL
		assertEquals(bps.get(1).bp, new BreakpointSummary(0, FWD, 4, 0, BWD, 7)); // INS
		assertEquals(bps.get(2).bp, new BreakpointSummary(0, FWD, 8, 0, BWD, 9)); // INV
		assertEquals(bps.get(3).bp, new BreakpointSummary(0, FWD, 10, 0, BWD, 11)); // INV
		assertEquals(bps.get(4).bp, new BreakpointSummary(0, FWD, 14, 0, BWD, 15)); // DUP
	}
	@Test
	public void should_output_variant_sequence() throws IOException {
		SimpleVariantChromosome svc = new SimpleVariantChromosome(getContext(), "polyACGT", 1, 0, true);
		svc.assemble(input, output, outputBedpe, false, Lists.newArrayList(SvType.DEL, SvType.INS, SvType.INV, SvType.DUP), Lists.newArrayList(2), 1, false);
		// 123456  78901234  5
		// ACGTAC  GTACGTAC  G
		// ||  ||  ||  ||    |
		// AC--ACNNGTGTGTACACG
		//  *|| *|| *|| *||||
		//   v1  v2  v3   v4
		List<String> fa = Files.readLines(input, StandardCharsets.US_ASCII);
		assertEquals("ACACCGGTGTGTACACG", fa.get(1));
	}
	@Test
	public void should_output_corresponding_reference_sequences_only() throws IOException {
		SimpleVariantChromosome svc = new SimpleVariantChromosome(getContext(), "polyACGT", 1, 0, true);
		svc.assemble(input, output, outputBedpe, true, Lists.newArrayList(SvType.DEL, SvType.INS, SvType.INV, SvType.DUP), Lists.newArrayList(2), 1, false);
		// 123456  78901234  5
		// ACGTAC  GTACGTAC  G
		// ||  ||  ||  ||    |
		// AC--ACNNGTGTGTACACG
		//  *|| *|| *|| *||||
		//   v1  v2  v3   v4
		List<String> fa = Files.readLines(input, StandardCharsets.US_ASCII);
		assertEquals(15, fa.get(3).length());
	}
	//@Test
	public void shouldCalculateMicrohomologyLength() throws IOException {
		// 12345678901234567890
		// NNANAAANAAAAAAANAAAAAAAAAAAAAAANAAAAAAAAAAAAAAAAAAAAAAAAAAA
		//           *|      *|
		//           v1      v2
		SimpleVariantChromosome svc = new SimpleVariantChromosome(getContext(), "Npower2", 2, 0, true);
		svc.assemble(input, output, outputBedpe, false, Lists.newArrayList(SvType.DEL), Lists.newArrayList(1), 5, false);
		List<IdsvVariantContext> vcf = super.getVcf(output, null);
		assertEquals("6", vcf.get(0).getAttribute("HOMLEN"));
		assertEquals(11, vcf.get(0).getStart());
		assertEquals(Lists.newArrayList("-3", "3"), vcf.get(0).getAttribute("CIPOS"));
		assertEquals(19, vcf.get(1).getStart());
		assertEquals("14", vcf.get(1).getAttribute("HOMLEN"));
		assertEquals(Lists.newArrayList("-3", "11"), vcf.get(1).getAttribute("CIPOS"));
	}
	@Test
	public void shouldFilterWhenVariantMatchesReferenceSequence() throws IOException {
		// 12345678
		// ACGTACGT
		//   *||
		SimpleVariantChromosome svc = new SimpleVariantChromosome(getContext(), "polyACGT", 0, 0, true);
		svc.assemble(input, output, outputBedpe, false, Lists.newArrayList(SvType.INV), Lists.newArrayList(2), 1, false);
		List<IdsvVariantContext> vcf = super.getVcf(output, null);
		assertEquals(1, vcf.get(0).getStart());
		assertTrue(vcf.get(0).getFilters().contains("REF"));
	}
}
