package au.edu.wehi.idsv.debruijn;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import com.google.common.collect.ImmutableList;

import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.debruijn.positional.KmerPathNode;


public class DeBruijnSequenceGraphNodeUtilTest extends TestHelper {
	@Test
	public void basesDifferent_should_count_bases_on_path() {
		BasePathGraph pg = PG(G(16)
				.add("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATC")
				.add("CATTAATCGCAAGAGCGGGTTGTATTCGACGTTAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATC")
				//                                   ^^
				);
		assertEquals(2, DeBruijnSequenceGraphNodeUtil.basesDifferent(pg.getK(), ImmutableList.of(pg.get("TCGACGTTAAGTCAGC")), ImmutableList.of(pg.get("TCGACGCCAAGTCAGC"))));
		
		
		String[] s = new String[] {  S(RANDOM).substring(0, 100), S(RANDOM).substring(200, 300) };
		pg = PG(G(20)
			.add(s[0])
			.add(s[1]));
		int mismatches = 0;
		for (int i = 0; i < 100; i++) if (s[0].charAt(i) != s[1].charAt(i)) mismatches++;
		assertEquals(mismatches, DeBruijnSequenceGraphNodeUtil.basesDifferent(pg.getK(),
				ImmutableList.of(pg.get(s[0].substring(0, 20))),
				ImmutableList.of(pg.get(s[1].substring(0, 20)))));
	}
	@Test
	public void split_bases_should_match_bases() {
		KmerPathNode a1 = KPN(25, "GATCGAGACCACGGTGAAACCCCGTTTCTATTAAAAATACAAAAAATTAGCCGGGCGCGGTGGCGGGCGCCTGTAGTCCCAG", 1, 1, false);
		KmerPathNode b1 = KPN(25, "GATCGAGACCACGGTGAAACCCCGTTTCTATTAAAA", 1, 1, false);
		KmerPathNode b2 = KPN(25,             "GGTGAAACCCCGTTTCTATTAAAAATACAAAAAATTAGCCGGGCGCGGTGGCGGGCGCCTGTAGTCCCAG", 13, 13, false);
		assertEquals(0, DeBruijnSequenceGraphNodeUtil.basesDifferent(25, ImmutableList.of(a1), ImmutableList.of(b1)));
		assertEquals(0, DeBruijnSequenceGraphNodeUtil.basesDifferent(25, ImmutableList.of(a1), ImmutableList.of(b1, b2)));
		assertEquals(0, DeBruijnSequenceGraphNodeUtil.reverseBasesDifferent(25, ImmutableList.of(a1), ImmutableList.of(b2)));
		assertEquals(0, DeBruijnSequenceGraphNodeUtil.reverseBasesDifferent(25, ImmutableList.of(a1), ImmutableList.of(b1, b2)));
	}
	/*
	@Test
	public void reverse_multisplit_bases_should_match_bases() {
		//                                                *                        *
		KmerPathNode a1 = KPN(25,                           "GTTTACAGCCTGTCTTATATCCCGAATAACGCACCGCCTATTCGAAC", 1, 1, false);
		KmerPathNode a2 = KPN(25, "GATTGACGTATCACAAGCCGGATGTTGTTTACAGCCTGTCTTATATCC", 1, 1, false);
		KmerPathNode b1 = KPN(25,                        "TTTGTTTACAGCCTGTCTTATATCCAGAATAACGCACCGCCTATTCGAAC", 1, 1, false);
		KmerPathNode b2 = KPN(25,                       "TTTTGTTTACAGCCTGTCTTATTCC", 1, 1, false);
		KmerPathNode b3 = KPN(25,                      "ATTTTGTTTACAGCCTGTCTTATTC", 1, 1, false);
		KmerPathNode b4 = KPN(25,                     "GATTTTGTTTACAGCCTGTCTTATT", 1, 1, false);
		KmerPathNode b5 = KPN(25,                    "GATTTTGTTTACAGCCTGTCTTATA", 1, 1, false);
		KmerPathNode b6 = KPN(25,                   "GGATTTTGTTTACAGCCTGTCTTAT", 1, 1, false);
		KmerPathNode b7 = KPN(25,                  "CGGATTTTGTTTACAGCCTGTCTTA", 1, 1, false);
		KmerPathNode b8 = KPN(25,                 "CCGGATTTTGTTTACAGCCTGTCTT", 1, 1, false);
		KmerPathNode b9 = KPN(25,                "GCCGGATTTTGTTTACAGCCTGTCT", 1, 1, false);
		KmerPathNode b10 = KPN(25,              "AGCCGGATTTTGTTTACAGCCTGTC", 1, 1, false);
		KmerPathNode b11 = KPN(25,             "AAGCCGGATTTTGTTTACAGCCTGT", 1, 1, false);
		KmerPathNode b12 = KPN(25,            "CAAGCCGGATTTTGTTTACAGCCTG", 1, 1, false);
		KmerPathNode b13 = KPN(25,           "ACAAGCCGGATTTTGTTTACAGCTT", 1, 1, false);
	}
	*/
}
