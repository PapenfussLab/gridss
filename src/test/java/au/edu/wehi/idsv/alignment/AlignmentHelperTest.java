package au.edu.wehi.idsv.alignment;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;



public class AlignmentHelperTest extends TestHelper {
	@Test
	public void alignment_should_ignore_case() {
		String ref = "gtacC";
		String seq = "GTACC";
		Alignment alignment = AlignerFactory.create().align_smith_waterman(B(seq), B(ref));
		assertEquals("5M", alignment.getCigar());
	}
	@Test
	public void test778_gridss_chr15_67104459_67104459f_108952510026_r() {
		String chr15_67104400_67104500 = "tgcagtttcttcctagcattgatggtctttacaatttggcatgtttttgcagtggctgggaccagttgttcctttccatgtttagtgcttccttcaggagc";
		String assembly = 					                               "TTTGGCATGTTTTTGCAGTGGCTGGGGGGGGGTGGTTTTT";
		Alignment alignment = AlignerFactory.create().align_smith_waterman(B(assembly), B(chr15_67104400_67104500));
		assertEquals("26M14S", alignment.getCigar());
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
		Alignment alignment = AlignerFactory.create().align_smith_waterman(B(assembly), B(chr9_gl000198_random_60000_60500));
		// Should match blastn result
		assertEquals("4S78M45S", alignment.getCigar());
	}
}
