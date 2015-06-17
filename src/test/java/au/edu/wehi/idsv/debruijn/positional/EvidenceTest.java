package au.edu.wehi.idsv.debruijn.positional;

import static org.junit.Assert.*;

import org.junit.Ignore;
import org.junit.Test;

import au.edu.wehi.idsv.SoftClipEvidence;
import au.edu.wehi.idsv.TestHelper;


public class EvidenceTest extends TestHelper {
	@Test
	public void softclip() {
		for (SoftClipEvidence sce : new SoftClipEvidence[] {
				SCE(FWD, withQual(new byte[] { 0,1,2,3,4,5,6,7,8,9,10}, withSequence("ACGTTATACCG", Read(0, 2, "1S4M6S")))),
				SCE(BWD, withQual(new byte[] { 0,1,2,3,4,5,6,7,8,9,10}, withSequence("ACGTTATACCG", Read(0, 2, "1S4M6S"))))}) {
			Evidence e = Evidence.create(4, sce, false);
			assertEquals(1, e.startPosition());
			assertEquals(1, e.endPosition());
			assertEquals(11-4+1, e.length());
			assertEquals(1, e.weight(0));
			assertEquals(2, e.weight(1));
			assertEquals(3, e.weight(2));
			assertEquals(4, e.weight(3));
			assertEquals(5, e.weight(4));
			assertEquals(6, e.weight(5));
			assertEquals(7, e.weight(6));
			assertEquals(8, e.weight(7));			
			assertEquals("ACGT", K(4, e.kmer(0)));
			assertEquals("CGTT", K(4, e.kmer(1)));
			assertEquals("GTTA", K(4, e.kmer(2)));
			assertEquals("TTAT", K(4, e.kmer(3)));
			assertEquals("TATA", K(4, e.kmer(4)));
			assertEquals("ATAC", K(4, e.kmer(5)));
			assertEquals("TACC", K(4, e.kmer(6)));
			assertEquals("ACCG", K(4, e.kmer(7)));
			assertFalse(e.isAnchored(0));
			assertTrue(e.isAnchored(1));
			assertFalse(e.isAnchored(2));
			assertFalse(e.isAnchored(3));
			assertFalse(e.isAnchored(4));
			assertFalse(e.isAnchored(5));
			assertFalse(e.isAnchored(6));
			assertFalse(e.isAnchored(7));
			assertEquals(sce.getEvidenceID(), e.evidenceId());
			for (int i = 0; i < e.length(); i++) {
				assertEquals(1 + i, e.node(i).startPosition());
				assertEquals(1 + i, e.node(i).endPosition());
				assertEquals(e.isAnchored(i), e.node(i).isReference());
				assertEquals(e.kmer(i), e.node(i).kmer());
				assertEquals(e.weight(i), e.node(i).weight());
				assertEquals(e, e.node(i).evidence());
			}
		}
	}
	@Test
	public void pair_fwd() {
		//          1         2         3
		// 123456789012345678901234567890123456789012345678901234567890
		// SMMMMMMMMMS                
		// |------------------| min
		//          MMMMMMMMMMM
		// |----------------------------| max
		//                    MMMMMMMMMMM
		MockSAMEvidenceSource ses = SES(20, 30);
		Evidence e = Evidence.create(4, NRRP(ses, withQual(new byte[] { 0,1,2,3,4,5,6,7,8,9,10}, withSequence("ACGTTATACCG", DP(0, 2, "1S9M1S", true, 1, 1, "11M", false)))));
		assertEquals(10, e.startPosition());
		assertEquals(20, e.endPosition());
		assertEquals(11-4+1, e.length());
		assertEquals("ACGT", K(4, e.kmer(0)));
		assertEquals("CGTT", K(4, e.kmer(1)));
		assertEquals("GTTA", K(4, e.kmer(2)));
		assertEquals("TTAT", K(4, e.kmer(3)));
		assertEquals("TATA", K(4, e.kmer(4)));
		assertEquals("ATAC", K(4, e.kmer(5)));
		assertEquals("TACC", K(4, e.kmer(6)));
		assertEquals("ACCG", K(4, e.kmer(7)));
		assertEquals(1, e.weight(0));
		assertEquals(2, e.weight(1));
		assertEquals(3, e.weight(2));
		assertEquals(4, e.weight(3));
		assertEquals(5, e.weight(4));
		assertEquals(6, e.weight(5));
		assertEquals(7, e.weight(6));
		assertEquals(8, e.weight(7));
		for (int i = 0; i < e.length(); i++) assertFalse(e.isAnchored(i));
	}
	@Test
	public void pair_bwd() {
		//          1         2         3
		// 123456789012345678901234567890123456789012345678901234567890
		//                    SMMMMMMMMMS                
		//           |------------------| min
		//           MMMMMMMMMMM
		// |----------------------------| max
		// MMMMMMMMMMM
		MockSAMEvidenceSource ses = SES(20, 30);
		Evidence e = Evidence.create(4, NRRP(ses, withQual(new byte[] { 0,1,2,3,4,5,6,7,8,9,10}, withSequence("ACGTTATACCG", DP(0, 21, "1S9M1S", false, 1, 1, "11M", true)))));
		assertEquals(1, e.startPosition());
		assertEquals(11, e.endPosition());
		assertEquals(11-4+1, e.length());
		assertEquals("ACGT", K(4, e.kmer(0)));
		assertEquals("CGTT", K(4, e.kmer(1)));
		assertEquals("GTTA", K(4, e.kmer(2)));
		assertEquals("TTAT", K(4, e.kmer(3)));
		assertEquals("TATA", K(4, e.kmer(4)));
		assertEquals("ATAC", K(4, e.kmer(5)));
		assertEquals("TACC", K(4, e.kmer(6)));
		assertEquals("ACCG", K(4, e.kmer(7)));
		assertEquals(1, e.weight(0));
		assertEquals(2, e.weight(1));
		assertEquals(3, e.weight(2));
		assertEquals(4, e.weight(3));
		assertEquals(5, e.weight(4));
		assertEquals(6, e.weight(5));
		assertEquals(7, e.weight(6));
		assertEquals(8, e.weight(7));
		for (int i = 0; i < e.length(); i++) assertFalse(e.isAnchored(i));
	}
	@Test
	public void should_reverse_comp_pair_if_required() {
		MockSAMEvidenceSource ses = SES(20, 30);
		Evidence e = Evidence.create(4, NRRP(ses, withQual(new byte[] { 0,1,2,3,4,5,6,7,8,9,10}, withSequence("ACGTTATACCG", DP(0, 21, "1S9M1S", false, 1, 1, "11M", false)))));
		assertEquals("CGGT", K(4, e.kmer(0))); // reverse comp of ending ACCG
	}
	@Test
	public void softclip_should_trim_soft_clip_on_other_side() {
		assertEquals(4, Evidence.create(2, SCE(FWD, Read(0, 1, "10H1S2M3S")), true).length());
		assertEquals(2, Evidence.create(2, SCE(BWD, Read(0, 1, "10H1S2M3S")), true).length());
	}
	@Test
	public void should_exclude_ambiguous_kmers() {
		Evidence e = Evidence.create(4, SCE(FWD, withQual(new byte[] { 0,1,2,3,4,5,6,7,8,9,10}, withSequence("ACGTNATACCG", Read(0, 2, "1S4M6S")))), false);
		assertNotNull(e.node(0));
		assertNull(e.node(1));
		assertNull(e.node(2));
		assertNull(e.node(3));
		assertNull(e.node(4));
		assertNotNull(e.node(5));
		e = Evidence.create(4, SCE(FWD, withQual(new byte[] { 0,1,2,3,4}, withSequence("ACGTN", Read(0, 2, "4M1S")))), false);
		assertNotNull(e.node(0));
		assertNull(e.node(1));
		e = Evidence.create(4, SCE(FWD, withQual(new byte[] { 0,1,2,3,4,5,6,7,8,9,10}, withSequence("NCGTTATACCN", Read(0, 2, "1S4M6S")))), false);
		assertNull(e.node(0));
		assertNotNull(e.node(1));
	}
	@Test
	@Ignore
	public void should_handle_cigar_indels() {
		// TODO: make sure both sides are at the correct position for cigars such as 5S1M1D1M5S
	}
}
