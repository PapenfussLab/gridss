package au.edu.wehi.idsv.sam;

import static org.junit.Assert.assertEquals;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;

import java.util.List;

import org.junit.Test;

import com.google.common.collect.Lists;


public class CigarUtilTest {
	public static List<CigarElement> C(String cigar) {
		return TextCigarCodec.getSingleton().decode(cigar).getCigarElements();
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
}
