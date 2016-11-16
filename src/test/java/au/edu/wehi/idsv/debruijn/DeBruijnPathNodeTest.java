package au.edu.wehi.idsv.debruijn;

import static org.junit.Assert.assertTrue;

import org.junit.Test;

import com.google.common.collect.ImmutableList;

import au.edu.wehi.idsv.TestHelper;


public class DeBruijnPathNodeTest extends TestHelper {
	@Test
	public void mergeShouldUpdateRefCounts() {
		DeBruijnGraphBase<DeBruijnNodeBase> g = new BaseGraph(3);
		DeBruijnNodeBase node1 = new DeBruijnNodeBase(0, 1, 0, "", true, 0, null);
		DeBruijnNodeBase node2 = new DeBruijnNodeBase(0, 1, 0, "", false, 0, null);
		DeBruijnPathNode<DeBruijnNodeBase> pn1 = new DeBruijnPathNode<DeBruijnNodeBase>(ImmutableList.of(node1), g);
		DeBruijnPathNode<DeBruijnNodeBase> pn2 = new DeBruijnPathNode<DeBruijnNodeBase>(ImmutableList.of(node2), g);
		assertTrue(pn1.isReference());
		pn1.merge(ImmutableList.of(pn2), g);
		assertTrue(pn1.isReference());
	}
}
