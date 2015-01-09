package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import jaligner.Alignment;
import jaligner.Sequence;

import org.junit.Test;


public class AlignmentHelperTest {
	@Test
	public void alignment_should_ignore_case() {
		String ref = "gtacC";
		String seq = "GTACC";
		Alignment alignment = AlignmentHelper.align_local(new Sequence(ref), new Sequence(seq));
		assertEquals("5M", AlignmentHelper.alignmentToCigar(alignment).toString());
	}
	@Test
	public void test778_gridss_chr15_67104459_67104459f_108952510026_r() {
		String chr15_67104400_67104500 = "tgcagtttcttcctagcattgatggtctttacaatttggcatgtttttgcagtggctgggaccagttgttcctttccatgtttagtgcttccttcaggagc";
		String assembly = 					                               "TTTGGCATGTTTTTGCAGTGGCTGGGGGGGGGTGGTTTTT";
		Alignment alignment = AlignmentHelper.align_local(new Sequence(chr15_67104400_67104500), new Sequence(assembly));
		assertEquals("26M14S", AlignmentHelper.alignmentToCigar(alignment).toString());
	}
	@Test
	public void test778_chr9_gl000198_random_60298() {
		String chr9_gl000198_random_60000_60500 =
			"cctccttaggatgcaaggttggtccagcatgtgcaaatcaataaatgtaatacatcacat" +
			"aaacagaactaaaaacaaaaaacagacgattatcccaatagatgcagaaaataataaaag" +
			"ccgtgtataacaaacccacagccaacaatctgttgaatgggcaaaagttgggagcgtttt" +
			"ttttttgaaaaccagcacaaggcaagaatgccctctctcaccatttctattcaacatagt" +
			"attgaaagtcctgaccagggcaatcaggctagagaaagaaatacaaaaagcatccaaata" +
			"ggaagagagaaagccaaactatttctgtttgcagataacataattctatatctagaaaac" +
			"cctgtagtttcagctacaaagcttctttagttaacaaacaacttcagcgaattttcaaga" +
			"tactaaatcaatgtgcaaagatcactaacatttctctacaccaatgataaccacactgac" +
			"agcaaaatcagaaaggccatc";
		String assembly = "TAGCAAGAGAGAAAGCCAAACTATTTCTGTTTGCAGATAACATAATTCTATATCTAGAAAACCCTGTAGTTTCAGGTACAAAACCCTGTAGTTTCAGGTACAATGTGCAAAGATCACTAACATTTCT";
		Alignment alignment = AlignmentHelper.align_local(new Sequence(chr9_gl000198_random_60000_60500), new Sequence(assembly));
		// Should match blastn result
		assertEquals("4S78M45S", AlignmentHelper.alignmentToCigar(alignment).toString());
	}
}
