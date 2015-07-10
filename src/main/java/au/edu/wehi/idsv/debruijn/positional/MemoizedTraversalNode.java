package au.edu.wehi.idsv.debruijn.positional;

public class MemoizedTraversalNode extends TraversalNode {
	private boolean valid = true;
	public boolean isValid() { return valid; }
	public void invalidate() {
		valid = false;
	}
	public MemoizedTraversalNode(KmerPathSubnode node, int baseScore) {
		super(node, baseScore);
	}
	public MemoizedTraversalNode(TraversalNode prev, KmerPathSubnode node) {
		super(prev, node);
	}
	public MemoizedTraversalNode(TraversalNode node, int start, int end) {
		super(node, start, end);
	}
}
