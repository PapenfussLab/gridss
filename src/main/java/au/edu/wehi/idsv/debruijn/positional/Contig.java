package au.edu.wehi.idsv.debruijn.positional;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Ordering;

import java.util.ArrayDeque;

public class Contig {
	public Contig(TraversalNode node, boolean hasReferenceSuccessor) {
		this.node = node;
		this.score = node.score + (hasReferenceSuccessor ? NonReferenceContigAssembler.ANCHORED_SCORE : 0);			
	}
	public Contig(TraversalNode node) {
		this.node = node;
		this.score = node.score;			
	}
	/**
	 * terminal contig node
	 */
	public final TraversalNode node;
	/**
	 * Final score for contig (including any anchor scoring bonus)
	 */
	public final int score;
	public ArrayDeque<KmerPathSubnode> toSubnodePath() {
		return node.toSubnodeNextPath();
	}
	@Override
	public String toString() {
		return String.format("Path Score %d, %s", score, node);
	}
	public static final Ordering<Contig> ByScoreDescPosition = new Ordering<Contig>() {
		@Override
		public int compare(Contig left, Contig right) {
			return ComparisonChain.start()
					.compare(right.score, left.score)
					.compare(left.node.node.lastStart(), right.node.node.lastStart())
					.result();
		}};
}