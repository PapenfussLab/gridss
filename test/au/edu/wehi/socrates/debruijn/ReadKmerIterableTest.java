package au.edu.wehi.socrates.debruijn;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

import java.util.Iterator;

import org.junit.Test;

import au.edu.wehi.socrates.TestHelper;

public class ReadKmerIterableTest extends TestHelper {
	@Test
	public void iterator_should_traverse_read_kmers() {
		ReadKmerIterable rki = new ReadKmerIterable(4, B("ACGTACGTT"), B("ABCDEAAAA"));
		Iterator<ReadKmer> i = rki.iterator();
		ReadKmer k;
		k = i.next(); assertEquals("ACGT", S(KmerEncodingHelper.encodedToPicardBases(k.kmer, 4)));
		k = i.next(); assertEquals("CGTA", S(KmerEncodingHelper.encodedToPicardBases(k.kmer, 4)));
		k = i.next(); assertEquals("GTAC", S(KmerEncodingHelper.encodedToPicardBases(k.kmer, 4)));
		k = i.next(); assertEquals("TACG", S(KmerEncodingHelper.encodedToPicardBases(k.kmer, 4)));
		k = i.next(); assertEquals("ACGT", S(KmerEncodingHelper.encodedToPicardBases(k.kmer, 4)));
		k = i.next(); assertEquals("CGTT", S(KmerEncodingHelper.encodedToPicardBases(k.kmer, 4)));
		assertFalse(i.hasNext());
	}
	@Test
	public void iterator_should_weight_kmer_by_min_base_quality_plus_one() {
		ReadKmerIterable rki = new ReadKmerIterable(4, B("ACGTACGTT"), B("AABCDEFGA"));
		Iterator<ReadKmer> i = rki.iterator();
		ReadKmer k;
		k = i.next(); assertEquals((byte)'A' + 1, k.weight);
		k = i.next(); assertEquals((byte)'A' + 1, k.weight);
		k = i.next(); assertEquals((byte)'B' + 1, k.weight);
		k = i.next(); assertEquals((byte)'C' + 1, k.weight);
		k = i.next(); assertEquals((byte)'D' + 1, k.weight);
		k = i.next(); assertEquals((byte)'A' + 1, k.weight);
		assertFalse(i.hasNext());
	}
}
