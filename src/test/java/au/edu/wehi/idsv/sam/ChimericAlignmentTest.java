package au.edu.wehi.idsv.sam;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.junit.Test;


public class ChimericAlignmentTest {
	@Test
	public void should_ignore_null_empty() {
		assertEquals(0, ChimericAlignment.getChimericAlignments((String)null).size());
		assertEquals(0, ChimericAlignment.getChimericAlignments("").size());
	}
	@Test
	public void should_decode_bwa_split_read() {
		String saz = "chr18,107870,-,8817S631M318S,30,39;chr18,108695,-,7874S237M1D203M1I5M2D215M1D40M1191S,0,48;chrUn_gl000216,155097,+,927S271M8568S,19,12;chrUn_gl000225,85612,-,7502S111M2153S,27,2;chrUn_gl000216,10548,-,724S76M8966S,0,1;chrY,13833846,+,7104S60M2602S,15,0;chrY,13869818,+,7747S58M1961S,0,0;chr2,89875994,+,4240S52M5474S,30,0;";
		List<ChimericAlignment> list = ChimericAlignment.getChimericAlignments(saz);
		assertEquals(8, list.size());
		assertEquals("chr18", list.get(0).rname);
		assertEquals(107870, list.get(0).pos);
		assertTrue(list.get(0).isNegativeStrand);
		assertEquals("8817S631M318S", list.get(0).cigar.toString());
		assertEquals(30, list.get(0).mapq);
		assertEquals(39, (int)list.get(0).nm);
	}
	@Test
	public void should_allow_missing_nm() {
		String sa = "chr18,107870,-,8817S631M318S,30,,";
		assertEquals(null, new ChimericAlignment(sa).nm);
	}
}
