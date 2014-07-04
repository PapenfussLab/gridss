package au.edu.wehi.idsv.debruijn.subgraph;

import static org.junit.Assert.*;

import org.junit.Test;


public class DeBruijnReadGraphTest {
	@Test
	public void reference_kmer_relevance_should_be_at_corresponding_reference_position() {
		fail("Currently not offsetting from anchor position - will call large backward reference subgraphs when they are still in scope");
		// workaround: we can pad out resolution by an additional read length - this will give us a closer anchor
	}
}
