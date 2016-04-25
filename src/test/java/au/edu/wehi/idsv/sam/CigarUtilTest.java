package au.edu.wehi.idsv.sam;

import static org.junit.Assert.assertEquals;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import org.apache.commons.lang3.tuple.Pair;
import org.junit.Test;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;


public class CigarUtilTest {
	public static List<CigarElement> C(String cigar) {
		return TextCigarCodec.decode(cigar).getCigarElements();
	}
	@Test
	public void encodeNegativeDeletion_should_convert_negD_to_xPxNxP() {
		assertEquals("1M5P5N5P2X1P1N1P", new Cigar(CigarUtil.encodeNegativeDeletion(Lists.newArrayList(
				new CigarElement(1, CigarOperator.M),
				new CigarElement(-5, CigarOperator.D),
				new CigarElement(2, CigarOperator.X),
				new CigarElement(-1, CigarOperator.D)))).toString());
	}
	@Test
	public void decodeNegativeDeletion_should_convert_xPxNxP_to_negD() {
		assertEquals("1M-5D2X-1D1M1P2N2P", new Cigar(CigarUtil.decodeNegativeDeletion(C("1M5P5N5P2X1P1N1P1M1P2N2P"))).toString());
	}
	@Test
	public void decodeNegativeDeletion_should_not_convert_PNP_if_mismatched_lengths() {
		assertEquals("1P2N2P", new Cigar(CigarUtil.decodeNegativeDeletion(C("1P2N2P"))).toString());
	}
	@Test
	public void readLength_should_match_read_length() {
		assertEquals(7, CigarUtil.readLength(C("1M2I5D2M2S1H")));
	}
	@Test
	public void referenceLength_should_match_reference_length() {
		assertEquals(5, CigarUtil.referenceLength(C("1M5I2D2X")));
	}
	private void assertSplitAfterReadPosition(String cigar, int position, String leftCigar, String rightCigar) {
		Pair<Cigar, Cigar> split = CigarUtil.splitAfterReadPosition(C(cigar), position);
		assertEquals(leftCigar, split.getLeft().toString());
		assertEquals(rightCigar, split.getRight().toString());
	}
	@Test
	public void splitAfterReadPosition_should_split_deletion() {
		assertSplitAfterReadPosition("1M1D1M", 0, "1M1S", "1S1M");
	}
	@Test
	public void splitAfterReadPosition_should_split_cigar() {
		assertSplitAfterReadPosition("4M", 0, "1M3S", "1S3M");
		assertSplitAfterReadPosition("4M", 1, "2M2S", "2S2M");
		assertSplitAfterReadPosition("4M", 2, "3M1S", "3S1M");
	}
	@Test
	public void splitAfterReadPosition_should_clean_cigar() {
		assertSplitAfterReadPosition("2M2D2M", 1, "2M2S", "2S2M");
		assertSplitAfterReadPosition("2M2I2M", 1, "2M4S", "4S2M");
		assertSplitAfterReadPosition("2M2I2M", 2, "2M4S", "4S2M");
		assertSplitAfterReadPosition("2M2I2M", 3, "2M4S", "4S2M");
	}
	@Test
	public void splitAfterReadPosition_should_decode_cigar() {
		assertSplitAfterReadPosition("2M10P10N10P2M", 1, "2M2S", "2S2M");
	}
	@Test
	public void offsetOf_should_decode_cigar() {
		assertEquals(-1, CigarUtil.offsetOf(new Cigar(C("2M3P3N3P2M")), 2));
	}
	@Test
	public void countMappedBases_should_count_MEX() {
		assertEquals(4 + 128 + 256, CigarUtil.countMappedBases(C("1H2H4M8I16D32N64P128=256X")));
	}
	@Test
	public void commonReferenceBases_should_count_overlap() {
		assertEquals(1, CigarUtil.commonReferenceBases(new Cigar(C("5M4S")), new Cigar(C("4S5M"))));
		assertEquals(7, CigarUtil.commonReferenceBases(new Cigar(C("1M2=4X")), new Cigar(C("1M2=4X"))));
	}
	@Test
	public void CigarOperatorIterator_should_match_cigar() {
		CigarOperatorIterator_test(C("1M"));
		CigarOperatorIterator_test(C("1M2M"));
		CigarOperatorIterator_test(C("1H2H4M8I16D32N64P128=256X"));
	}
	private void CigarOperatorIterator_test(List<CigarElement> l) {
		List<CigarOperator> co = new ArrayList<CigarOperator>();
		for (CigarElement e : l) {
			for (int i = 0; i < e.getLength(); i++) {
				co.add(e.getOperator());
			}
		}
		List<CigarOperator> result = Lists.newArrayList(new CigarUtil.CigarOperatorIterator(l));
		assertEquals(co, result);
	}
	@Test
	public void clean_should_merge_adjacent_and_trim_zero_size_operators() {
		assertEquals("4M", new Cigar(CigarUtil.clean(C("1M0I0D1M0M2M0S"))).toString());
	}
	@Test
	public void clean_should_fix_open_indels() {
		assertEquals("2S1M", new Cigar(CigarUtil.clean(C("1S1D1I1M"))).toString());
	}
	@Test
	public void offsetOf_should_use_first_alignment_at_or_after_offset_position() {
		assertEquals(0, CigarUtil.offsetOf(new Cigar(C("3S3M")), 0));
		assertEquals(0, CigarUtil.offsetOf(new Cigar(C("3S3M")), 1));
		assertEquals(0, CigarUtil.offsetOf(new Cigar(C("3S3M")), 2));
		assertEquals(0, CigarUtil.offsetOf(new Cigar(C("3S3M")), 3));
		assertEquals(1, CigarUtil.offsetOf(new Cigar(C("3S3M")), 4));
		assertEquals(2, CigarUtil.offsetOf(new Cigar(C("3S3M")), 5));
		assertEquals(3, CigarUtil.offsetOf(new Cigar(C("3S3M1I1M")), 6));
		assertEquals(3, CigarUtil.offsetOf(new Cigar(C("3S3M1I1M")), 7));
		assertEquals(6, CigarUtil.offsetOf(new Cigar(C("3S3M1I1M2D1M")), 8));
	}
	@Test
	public void offsetOf_should_respect_indels() {
		assertEquals(2, CigarUtil.offsetOf(new Cigar(C("1M1D1M")), 1));
	}
	@Test
	public void trimReadBases_should_remove_bases() {
		assertEquals("1S2M", CigarUtil.trimReadBases(new Cigar(C("3S3M")), 2, 1).toString());
		assertEquals("3M", CigarUtil.trimReadBases(new Cigar(C("4S1M2D3M")), 5, 0).toString());
		assertEquals("2S1M", CigarUtil.trimReadBases(new Cigar(C("4S1M2D3M")), 2, 3).toString());
	}
	@Test
	public void commonReferenceBases_should_count_common() {
		Cigar c1 = new Cigar(C("7S65M"));
		Cigar c2 = new Cigar(C("68M4S"));
		assertEquals(61, CigarUtil.commonReferenceBases(c1, c2));
		assertEquals(65, CigarUtil.countMappedBases(c1.getCigarElements()));
		assertEquals(68, CigarUtil.countMappedBases(c2.getCigarElements()));
	}
	@Test
	public void commonReferenceBases_should_match_unencoded_and_decoded_negative_deletions() {
		Cigar c1 = new Cigar(ImmutableList.of(
				new CigarElement(2, CigarOperator.M),
				new CigarElement(-2, CigarOperator.D),
				new CigarElement(2, CigarOperator.M)
				));
		assertEquals(4, CigarUtil.commonReferenceBases(c1, new Cigar(C("2M5P5N5P2M"))));
	}
	//@Test
	public void test_realign_anchor_direction_changing_na12878() throws IOException {
		Files.write(new File("C:/temp/asm.csv").toPath(), Files.lines(new File("C:/dev/idsv/src/test/resources/asm.csv").toPath())
			.map(line -> {
			String[] split = line.split(",");
			Cigar before = new Cigar(C(split[0]));
			int bases = CigarUtil.commonReferenceBases(before, new Cigar(C(split[1])));
			return line + "," + Float.toString(bases / (float)CigarUtil.countMappedBases(before.getCigarElements()));
		}).collect(Collectors.toList()));
		//dt <- read.csv("c:/temp/asm.csv", header=FALSE)
		//ggplot(dt) + aes(x=dt$V3) + geom_histogram()
	}
	@Test
	public void asUngapped_should_convert_operators() {
		assertEquals("10M", new Cigar(CigarUtil.asUngapped(new Cigar(C("10M")), true)).toString());
		assertEquals("10M5S", new Cigar(CigarUtil.asUngapped(new Cigar(C("10M5S")), true)).toString());
		assertEquals("1H2S10M5S", new Cigar(CigarUtil.asUngapped(new Cigar(C("1H2S10M5S")), true)).toString());
		assertEquals("10M15S", new Cigar(CigarUtil.asUngapped(new Cigar(C("10M1D10M5S")), true)).toString());
		assertEquals("10M16S", new Cigar(CigarUtil.asUngapped(new Cigar(C("10M1I10M5S")), true)).toString());
		assertEquals("11S10M5S", new Cigar(CigarUtil.asUngapped(new Cigar(C("10M1I10M5S")), false)).toString());
	}
}
