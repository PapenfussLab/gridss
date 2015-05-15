package au.edu.wehi.idsv.debruijn;

import java.util.List;

import org.junit.Test;
import static org.junit.Assert.assertEquals;

import com.google.common.collect.Lists;

import au.edu.wehi.idsv.TestHelper;


public class PackedReadKmerListTest extends TestHelper {
	public void matches(int k, String sequence, String weights) {
		matches(k, sequence, weights, false, false);
		matches(k, sequence, weights, false, true);
		matches(k, sequence, weights, true, false);
		matches(k, sequence, weights, true, true);
	}
	public void matches(int k, String sequence, String weights, boolean reverse, boolean complement) {
		ReadKmerIterable it = new ReadKmerIterable(k, B(sequence), B(weights), reverse, complement);
		PackedReadKmerList packed = new PackedReadKmerList(k, B(sequence), B(weights), reverse, complement);
		List<ReadKmer> rk = Lists.newArrayList(it.iterator());
		assertEquals(rk.size(), packed.length());
		for (int i = 0; i < rk.size(); i++) {
			assertEquals(rk.get(i).kmer, packed.kmer(i));
			assertEquals(rk.get(i).weight, packed.weight(i));
		}
	}
	@Test
	public void shouldMatchReadKmerIterable() {
		matches(1, "A", "1");
		matches(1, "AG", "12");
		matches(1, "AGT", "123");
		matches(1, "AGTC", "1234");
		matches(1, "AGTCA", "12345");
		matches(1, "AGTCAG", "123456");
		
		matches(2, "AG", "12");
		matches(2, "AGT", "123");
		matches(2, "AGTC", "1234");
		matches(2, "AGTCA", "12345");
		matches(2, "AGTCAG", "123456");
		
		matches(3, "AGT", "123");
		matches(3, "AGTC", "1234");
		matches(3, "AGTCA", "12345");
		matches(3, "AGTCAG", "123456");
		
		matches(4, "AGTC", "1234");
		matches(4, "AGTCA", "12345");
		matches(4, "AGTCAG", "123456");
		
		matches(8, "GTACGTAC", "12341234");
		matches(9, "GTACTGTACGTAC", "1234512341234");
		matches(10, "GTACTGTACGTAC", "1234512341234");
		matches(11, "GTACTGTACGTAC", "1234512341234");
		matches(12, "GTACTGTACGTAC", "1234512341234");
	}
	@Test
	public void shouldLongPack2bitEncodedBases() {
		matches(8, "AAAAAAAA", "AAAAAAAA", false, false);
	}
	@Test
	public void should_complement() {
		assertEquals("AGT", K(3, new PackedReadKmerList(3, B("TCA"), B("AAAA"), false, true).kmer(0)));
	}
	@Test
	public void should_reverse() {
		assertEquals("ACT", K(3, new PackedReadKmerList(3, B("TCA"), B("AAAA"), true, false).kmer(0)));
	}
	@Test
	public void should_reverse_complement() {
		assertEquals("TGA", K(3, new PackedReadKmerList(3, B("TCA"), B("AAA"), true, true).kmer(0)));
	}
}
