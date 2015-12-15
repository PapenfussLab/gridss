package au.edu.wehi.idsv.sam;

import static org.junit.Assert.assertEquals;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

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
	@Test
	public void splitAtLargestIndel_should_split_at_indel() {
		List<List<CigarElement>> split = CigarUtil.splitAtLargestIndel(C("1M2I1M"));
		assertEquals("1M", new Cigar(split.get(0)).toString());
		assertEquals("2I", new Cigar(split.get(1)).toString());
		assertEquals("1M", new Cigar(split.get(2)).toString());
	}
	@Test
	public void splitAtLargestIndel_should_include_successive_indel_ops_in_split() {
		List<List<CigarElement>> split = CigarUtil.splitAtLargestIndel(C("1M2I2D1M"));
		assertEquals("1M", new Cigar(split.get(0)).toString());
		assertEquals("2I2D", new Cigar(split.get(1)).toString());
		assertEquals("1M", new Cigar(split.get(2)).toString());
	}
	@Test
	public void splitAtLargestIndel_should_sort_on_largest_abs_event_size() {
		List<List<CigarElement>> split = CigarUtil.splitAtLargestIndel(C("1M5P5N5P2M3I"));
		assertEquals("1M", new Cigar(split.get(0)).toString());
		assertEquals("-5D", new Cigar(split.get(1)).toString());
		assertEquals("2M3I", new Cigar(split.get(2)).toString());
	}
	@Test
	public void splitAtLargestIndel_should_create_placeholder_0M_if_no_indel() {
		List<List<CigarElement>> split = CigarUtil.splitAtLargestIndel(C("10M"));
		assertEquals("10M", new Cigar(split.get(0)).toString());
		assertEquals("0M", new Cigar(split.get(2)).toString());
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
}
