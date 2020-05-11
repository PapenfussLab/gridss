package au.edu.wehi.idsv.sim;

import au.edu.wehi.idsv.sim.SequentialVariantPlacer.ContigExhaustedException;
import org.junit.Test;

import static org.junit.Assert.assertEquals;


public class SequentialVariantPlacerTest extends au.edu.wehi.idsv.TestHelper {
	@Test
	public void should_pad_contig_start() throws ContigExhaustedException {
		SequentialVariantPlacer svp = new SequentialVariantPlacer(POLY_A);
		svp.setDistanceBetweenVariants(1);
		svp.setMinimumDistanceFromN(1);
		assertEquals(2, svp.getNext(1));
		
		svp = new SequentialVariantPlacer(POLY_A);
		svp.setDistanceBetweenVariants(2);
		svp.setMinimumDistanceFromN(2);
		assertEquals(3, svp.getNext(1));
	}
	@Test(expected=ContigExhaustedException.class)
	public void should_pad_contig_end() throws ContigExhaustedException {
		SequentialVariantPlacer svp = new SequentialVariantPlacer(B("1234567890"));
		svp.setDistanceBetweenVariants(3);
		svp.setMinimumDistanceFromN(3);
		// 1234567890
		//    ^   ^ 
		svp.getNext(1);
		svp.getNext(1); // not enough padding to end of contig
	}
	@Test(expected=IllegalArgumentException.class)
	public void shouldRequirePositiveFeatureSize() throws ContigExhaustedException  {
		SequentialVariantPlacer svp = new SequentialVariantPlacer(POLY_A);
		svp.getNext(0);
	}
	@Test
	public void should_place_apart_according_to_padding() throws ContigExhaustedException  {
		SequentialVariantPlacer svp = new SequentialVariantPlacer(POLY_A);
		svp.setDistanceBetweenVariants(0);
		assertEquals(svp.getNext(1) + 1, svp.getNext(1)); // adjacent
		svp.setDistanceBetweenVariants(1);
		assertEquals(svp.getNext(1) + 2, svp.getNext(1)); // 1 base padding
		svp.setDistanceBetweenVariants(2);
		assertEquals(svp.getNext(1) + 3, svp.getNext(1)); // 2 base padding
	}
	@Test
	public void should_pad_from_end_of_last_variant() throws ContigExhaustedException  {
		SequentialVariantPlacer svp = new SequentialVariantPlacer(POLY_A);
		svp.setDistanceBetweenVariants(0);
		int pos = svp.getNext(1);
		assertEquals(pos + 1, pos = svp.getNext(2));
		assertEquals(pos + 2, pos = svp.getNext(3));
		assertEquals(pos + 3, pos = svp.getNext(4));
		assertEquals(pos + 4, pos = svp.getNext(1));
		
		svp.setDistanceBetweenVariants(1);
		assertEquals(pos + 2, pos = svp.getNext(2));
		assertEquals(pos + 3, pos = svp.getNext(3));
		assertEquals(pos + 4, pos = svp.getNext(4));
		assertEquals(pos + 5, pos = svp.getNext(1));
	}
	@Test
	public void should_N_pad_at_start_of_variant() throws ContigExhaustedException  {
		SequentialVariantPlacer svp = new SequentialVariantPlacer(B("NNNNNNAAAAAAAAAA"));
		//															 1234567890123456
		//															 ------21012-----
		svp.setDistanceBetweenVariants(0);
		svp.setMinimumDistanceFromN(2);
		assertEquals(9, svp.getNext(1));
	}
	@Test
	public void should_N_pad_at_end_of_variant() throws ContigExhaustedException  {
		SequentialVariantPlacer svp = new SequentialVariantPlacer(
			B("NAACCANAACCAANAAAAAAAAAAAAAAAAAAAAAAAA"));
		//	   1234567890123456
		svp.setDistanceBetweenVariants(0);
		svp.setMinimumDistanceFromN(2);
		assertEquals(10, svp.getNext(2)); // position 4 does not have enough N padding after the variant
	}
	@Test
	public void zero_N_pad_should_allow_variants_at_N() throws ContigExhaustedException  {
		SequentialVariantPlacer svp = new SequentialVariantPlacer(B("NNNNNNAAAAAAAAAA"));
		//															 1234567890123456
		//															 ------21012-----
		svp.setDistanceBetweenVariants(0);
		svp.setMinimumDistanceFromN(0);
		assertEquals(1, svp.getNext(1));
	}
	@Test
	public void zero_pad_should_place_variants_adjacently() throws ContigExhaustedException  {
		SequentialVariantPlacer svp = new SequentialVariantPlacer(POLY_A);
		svp.setDistanceBetweenVariants(0);
		int p1 = svp.getNext(1);
		int p2 = svp.getNext(1);
		assertEquals(p1 + 1, p2);
	}
	@Test
	public void minimumDistanceFromN_should_default_to_distanceBetweenVariants() {
		SequentialVariantPlacer svp = new SequentialVariantPlacer(POLY_A);
		assertEquals(svp.getDistanceBetweenVariants(), svp.getMinimumDistanceFromN());
		svp.setDistanceBetweenVariants(0);
		assertEquals(0, svp.getMinimumDistanceFromN());
		svp.setDistanceBetweenVariants(1);
		assertEquals(1, svp.getMinimumDistanceFromN());
	}
}
