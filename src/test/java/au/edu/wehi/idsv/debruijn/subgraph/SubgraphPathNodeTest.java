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
		DeBruijnReadGraph g = RG(4);
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
		DeBruijnReadGraph g = RG(4);
		g.addEvidence(SCE(FWD, withSequence("TTTTAAAAC", Read(0, 10, "4M5S"))));
		assertFalse(new SubgraphPathNode(ImmutableList.of(KmerEncodingHelper.picardBaseToEncoded(4, B("TTTT"))), g).containsNonReferenceKmer());
		assertTrue(new SubgraphPathNode(ImmutableList.of(KmerEncodingHelper.picardBaseToEncoded(4, B("TTTA"))), g).containsNonReferenceKmer()); 
		assertTrue(new SubgraphPathNode(ImmutableList.of(
				KmerEncodingHelper.picardBaseToEncoded(4, B("TTTT")),
				KmerEncodingHelper.picardBaseToEncoded(4, B("TTTA"))
				), g).containsNonReferenceKmer());  
	}
	@Test
	public void onKmersChanged_should_recalc_reference() {
		DeBruijnReadGraph g = RG(4);
		g.addEvidence(SCE(FWD, withSequence("TTTTAAAA", Read(0, 10, "4M4S"))));
		SubgraphPathNode r = new SubgraphPathNode(ImmutableList.of(KmerEncodingHelper.picardBaseToEncoded(4, B("TTTT"))), g);
		SubgraphPathNode n = new SubgraphPathNode(ImmutableList.of(KmerEncodingHelper.picardBaseToEncoded(4, B("AAAA"))), g);
		assertTrue(r.containsReferenceKmer());
		assertFalse(n.containsReferenceKmer());
		n.merge(ImmutableList.of(r), g);		
		assertTrue(n.containsReferenceKmer());
	}
}
