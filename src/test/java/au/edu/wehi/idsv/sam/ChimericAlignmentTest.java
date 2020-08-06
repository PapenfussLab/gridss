package au.edu.wehi.idsv.sam;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;
import htsjdk.samtools.TextCigarCodec;
import org.junit.Assert;
import org.junit.Test;

import java.util.List;

import static org.junit.Assert.*;


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
	@Test
	public void getAlignedBaseReadOffset_should_consider_strand() {
		// 0123456789
		// sMMssshhhh
		assertEquals(1, new ChimericAlignment(null, 0, false, TextCigarCodec.decode("1S2M3S4H"), 0, 0).getFirstAlignedBaseReadOffset());
		assertEquals(2, new ChimericAlignment(null, 0, false, TextCigarCodec.decode("1S2M3S4H"), 0, 0).getLastAlignedBaseReadOffset());
		// 0123456789
		// hhhhsssMMs
		assertEquals(7, new ChimericAlignment(null, 0, true, TextCigarCodec.decode("1S2M3S4H"), 0, 0).getFirstAlignedBaseReadOffset());
		assertEquals(8, new ChimericAlignment(null, 0, true, TextCigarCodec.decode("1S2M3S4H"), 0, 0).getLastAlignedBaseReadOffset());
	}
	@Test
	public void getAlignedBaseReadOffset_should_include_hard_clipping() {
		// 01234567
		//  MM
		assertEquals(1, new ChimericAlignment(null, 0, false, TextCigarCodec.decode("1H2M3H"), 0, 0).getFirstAlignedBaseReadOffset());
		assertEquals(2, new ChimericAlignment(null, 0, false, TextCigarCodec.decode("1H2M3H"), 0, 0).getLastAlignedBaseReadOffset());
		// 012345
		// hhhMMh
		assertEquals(3, new ChimericAlignment(null, 0, true, TextCigarCodec.decode("1H2M3H"), 0, 0).getFirstAlignedBaseReadOffset());
		assertEquals(4, new ChimericAlignment(null, 0, true, TextCigarCodec.decode("1H2M3H"), 0, 0).getLastAlignedBaseReadOffset());
	}
	@Test
	public void should_support_BELAN_format_with_HLA_types() {
		//chr:start|strand|cigar|mapq
		String bealn = "SN:HLA-DRB1*15:03:01:01:1-1:107870|-|8817S631M318S|30";
		ChimericAlignment ca = ChimericAlignment.parseBEALNAlignment(bealn);
		assertEquals("SN:HLA-DRB1*15:03:01:01:1-1", ca.rname);
		assertEquals(107870, ca.pos);
		assertTrue(ca.isNegativeStrand);
		assertEquals("8817S631M318S", ca.cigar.toString());
		assertEquals(30, ca.mapq);
		assertNull(ca.nm);
	}

	/**
	 * No mapq reported for XA tags
	 */
	@Test
	public void should_support_BELAN_XA() {
		String bealn = "SN:HLA-DRB1*15:03:01:01:1-1:107870|-|8817S631M318S|";
		ChimericAlignment ca = ChimericAlignment.parseBEALNAlignment(bealn);
		assertEquals(255, ca.mapq);
	}

	@Test
	public void intervals_should_merge_overlapping() {
		String saz = "polyA,1,+,5S10M,0,0;polyA,1,+,10M5S,0,0";
		List<ChimericAlignment> list = ChimericAlignment.getChimericAlignments(saz);
		RangeSet<Integer> rs = ChimericAlignment.getAlignedIntervals(list);
		Assert.assertEquals(TreeRangeSet.create(ImmutableList.of(
				Range.closedOpen(0, 15))),
				ChimericAlignment.getAlignedIntervals(list));
		Assert.assertEquals(TreeRangeSet.create(ImmutableList.of(
				)),
				ChimericAlignment.getUnalignedIntervals(list));
	}
	@Test
	public void intervals_should_merge_adjacent() {
		String saz = "polyA,1,+,7S8M,0,0;polyA,1,+,7M8S,0,0";
		List<ChimericAlignment> list = ChimericAlignment.getChimericAlignments(saz);
		RangeSet<Integer> rs = ChimericAlignment.getAlignedIntervals(list);
		Assert.assertEquals(TreeRangeSet.create(ImmutableList.of(
				Range.closedOpen(0, 15))),
				ChimericAlignment.getAlignedIntervals(list));
		Assert.assertEquals(TreeRangeSet.create(ImmutableList.of(
				)),
				ChimericAlignment.getUnalignedIntervals(list));
	}

	@Test
	public void intervals_should_union_aligned() {
		// 01234567
		// M MM MMM
		String saz = "polyA,1,+,1M7S,0,0;polyA,1,+,2S2M4S,0,0;polyA,1,+,5S3M,0,0;";
		List<ChimericAlignment> list = ChimericAlignment.getChimericAlignments(saz);
		RangeSet<Integer> rs = ChimericAlignment.getAlignedIntervals(list);
		Assert.assertEquals(TreeRangeSet.create(ImmutableList.of(
				Range.closedOpen(0, 1),
				Range.closedOpen(2, 4),
				Range.closedOpen(5, 8))),
				ChimericAlignment.getAlignedIntervals(list));
		Assert.assertEquals(TreeRangeSet.create(ImmutableList.of(
				Range.closedOpen(1, 2),
				Range.closedOpen(4, 5))),
				ChimericAlignment.getUnalignedIntervals(list));
	}
	@Test
	public void intervals_should_consider_strand() {
		// 01234567
		// >>
		//   <<<<<
		String saz = "polyA,1,+,2M6S,0,0;polyA,1,-,1S5M2S,0,0";
		List<ChimericAlignment> list = ChimericAlignment.getChimericAlignments(saz);
		RangeSet<Integer> rs = ChimericAlignment.getAlignedIntervals(list);
		Assert.assertEquals(TreeRangeSet.create(ImmutableList.of(
				Range.closedOpen(0, 7))),
				ChimericAlignment.getAlignedIntervals(list));
		Assert.assertEquals(TreeRangeSet.create(ImmutableList.of(
				Range.closedOpen(7, 8))),
				ChimericAlignment.getUnalignedIntervals(list));
	}

	@Test
	public void intervals_should_convert_hard_clipping() {
		String saz = "polyA,1,+,5H10M,0,0";
		List<ChimericAlignment> list = ChimericAlignment.getChimericAlignments(saz);
		RangeSet<Integer> rs = ChimericAlignment.getAlignedIntervals(list);
		Assert.assertEquals(TreeRangeSet.create(ImmutableList.of(
				Range.closedOpen(5, 15))),
				ChimericAlignment.getAlignedIntervals(list));
		Assert.assertEquals(TreeRangeSet.create(ImmutableList.of(
				Range.closedOpen(0, 5))),
				ChimericAlignment.getUnalignedIntervals(list));
	}
}
