package au.edu.wehi.idsv.debruijn;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Iterator;

import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;

public class ReadKmerIterableTest extends TestHelper {
	@Test
	public void iterator_should_traverse_read_kmers() {
		ReadKmerIterable rki = new ReadKmerIterable(4, B("ACGTACGTT"), B("ABCDEAAAA"));
		Iterator<ReadKmer> i = rki.iterator();
		ReadKmer k;
		k = i.next(); assertEquals("ACGT", S(KmerEncodingHelper.encodedToPicardBases(4, k.kmer)));
		k = i.next(); assertEquals("CGTA", S(KmerEncodingHelper.encodedToPicardBases(4, k.kmer)));
		k = i.next(); assertEquals("GTAC", S(KmerEncodingHelper.encodedToPicardBases(4, k.kmer)));
		k = i.next(); assertEquals("TACG", S(KmerEncodingHelper.encodedToPicardBases(4, k.kmer)));
		k = i.next(); assertEquals("ACGT", S(KmerEncodingHelper.encodedToPicardBases(4, k.kmer)));
		k = i.next(); assertEquals("CGTT", S(KmerEncodingHelper.encodedToPicardBases(4, k.kmer)));
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
	@Test
	public void iterator_should_track_ambiguous_bases() {
		ReadKmerIterable rki = new ReadKmerIterable(4, B("ACTTNNCAGTN"), B("AAAAAAAAAAA"));
		Iterator<ReadKmer> i = rki.iterator();
		ReadKmer k;
		k = i.next(); assertFalse(k.containsAmbiguousBases);
		k = i.next(); assertTrue(k.containsAmbiguousBases);
		k = i.next(); assertTrue(k.containsAmbiguousBases);
		k = i.next(); assertTrue(k.containsAmbiguousBases);
		k = i.next(); assertTrue(k.containsAmbiguousBases);
		k = i.next(); assertTrue(k.containsAmbiguousBases);
		k = i.next(); assertFalse(k.containsAmbiguousBases);
		k = i.next(); assertTrue(k.containsAmbiguousBases);
		assertFalse(i.hasNext());
	}
	@Test
	public void iterator_should_track_ambiguous_bases_in_start_kmer() {
		ReadKmerIterable rki = new ReadKmerIterable(4, B("ANAAAAT"), B("AAAAAAA"));
		Iterator<ReadKmer> i = rki.iterator();
		ReadKmer k;
		k = i.next(); assertTrue(k.containsAmbiguousBases);
		k = i.next(); assertTrue(k.containsAmbiguousBases);
		k = i.next(); assertFalse(k.containsAmbiguousBases);
		k = i.next(); assertFalse(k.containsAmbiguousBases);
		assertFalse(i.hasNext());
		
		i = new ReadKmerIterable(4, B("ANNAAAT"), B("AAAAAAA")).iterator();
		k = i.next(); assertTrue(k.containsAmbiguousBases);
		k = i.next(); assertTrue(k.containsAmbiguousBases);
		k = i.next(); assertTrue(k.containsAmbiguousBases);
		k = i.next(); assertFalse(k.containsAmbiguousBases);
		assertFalse(i.hasNext());
		
		i = new ReadKmerIterable(4, B("NAANAATA"), B("AAAAAAAA")).iterator();
		k = i.next(); assertTrue(k.containsAmbiguousBases);
		k = i.next(); assertTrue(k.containsAmbiguousBases);
		k = i.next(); assertTrue(k.containsAmbiguousBases);
		k = i.next(); assertTrue(k.containsAmbiguousBases);
		k = i.next(); assertFalse(k.containsAmbiguousBases);
		assertFalse(i.hasNext());
	}
}
