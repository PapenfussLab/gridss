package au.edu.wehi.idsv.debruijn;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;
import com.google.common.collect.ImmutableList;


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
}
