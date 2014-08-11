package au.edu.wehi.idsv.debruijn.subgraph;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;

import com.google.common.collect.ImmutableList;


public class SubgraphPathNodeTest extends TestHelper {
	@Test
	public void containsReferenceKmer() {
		// no ref kmers if anchor too short
		DeBruijnReadGraph g = G(4, FWD);
		g.addEvidence(SCE(FWD, withSequence("TTTTAAAAC", Read(0, 10, "4M5S"))));
		assertTrue(new SubgraphPathNode(ImmutableList.of(KmerEncodingHelper.picardBaseToEncoded(4, B("TTTT"))), g).containsReferenceKmer());
		assertFalse(new SubgraphPathNode(ImmutableList.of(KmerEncodingHelper.picardBaseToEncoded(4, B("TTTA"))), g).containsReferenceKmer()); 
		assertTrue(new SubgraphPathNode(ImmutableList.of(
				KmerEncodingHelper.picardBaseToEncoded(4, B("TTTT")),
				KmerEncodingHelper.picardBaseToEncoded(4, B("TTTA"))
				), g).containsReferenceKmer());  
	}
	@Test
	public void containsNonReferenceKmer() {
		// no ref kmers if anchor too short
		DeBruijnReadGraph g = G(4, FWD);
		g.addEvidence(SCE(FWD, withSequence("TTTTAAAAC", Read(0, 10, "4M5S"))));
		assertFalse(new SubgraphPathNode(ImmutableList.of(KmerEncodingHelper.picardBaseToEncoded(4, B("TTTT"))), g).containsNonReferenceKmer());
		assertTrue(new SubgraphPathNode(ImmutableList.of(KmerEncodingHelper.picardBaseToEncoded(4, B("TTTA"))), g).containsNonReferenceKmer()); 
		assertTrue(new SubgraphPathNode(ImmutableList.of(
				KmerEncodingHelper.picardBaseToEncoded(4, B("TTTT")),
				KmerEncodingHelper.picardBaseToEncoded(4, B("TTTA"))
				), g).containsNonReferenceKmer());  
	}
}