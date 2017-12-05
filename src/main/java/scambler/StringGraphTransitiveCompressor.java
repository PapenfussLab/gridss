package scambler;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

import au.edu.wehi.idsv.Defaults;
import scambler.SgNode.State;
import scambler.SgNode.TransitiveReductionMark;

/**
 * Streaming implementation of Myer 2005 String Graph transitive reduction algorithm
 * that also performs node compression.
 * 
 * Nodes are output when we can guarantee that no edges in or out of the node will
 * be modified after emission.
 * 
 * Transitive reduction can remove an edge coming into a node when any parent is processed
 * 
 * Transitive reduction can remove an edge coming out a node when the node itself is processed
 * 
 * Therefore, to emit, both the node and all parents must have had transitive reduction
 * performed on them.
 * 
 * To perform transitive reduction, both the node, and all children must have all
 * outgoing edges loaded.
 * 
 * Node compression: nodes can become unbranching either when
 * a) The load is loaded. This is the simple case and the node can immediate be compressed
 * b) Transitive reduction occurs. The node to compress can be either
 * bi) The node being reduced
 * bii) The child of the node being reduced
 * 
 * Therefore, we have the following constraints:
 * 
 * Transitive reduction:
 * u parent: >= loaded
 * v node: >= loaded
 * w child: >= loaded
 * x grandchild: no constraints (vw and wx exist since v,w loaded) 
 * 
 * 
 * @author Daniel Cameron
 *
 */
public class StringGraphTransitiveCompressor implements Iterator<SgNode> {
	private final Iterator<SgNode> it;
	//private final LinearGenomicCoordinate lgc;
	private final Deque<SgNode> output = new ArrayDeque<>();
	//private final Set<SgNode> canProcess = new HashSet<>();
	private final int fuzz;
	private final boolean compressUnbranchingPaths;
	public StringGraphTransitiveCompressor(Iterator<SgNode> it, int fuzz, boolean compressUnbranchingPaths) {
		this.it = it;
		this.fuzz = fuzz;
		this.compressUnbranchingPaths = compressUnbranchingPaths; 
	}
	private void ensureNext() {
		while (it.hasNext() && output.isEmpty()) {
			load(it.next());
		}
	}
	/**
	 * Transitive reduction requires that all children be fully constructed
	 * @param node
	 * @return
	 */
	private boolean isBlockingParentTransitiveReduction(SgNode node) {
		switch (node.state) {
		case Removed:
		case UnderConstruction:
			assert(node.state != State.Removed);
			return true;
		case AwaitingReduction:
		case Reduced:
		case EmittedFromReduction:
		default:
			return false;
		}
	}
	/**
	 * Emission can occur when no edges involving that node can change
	 * 
	 * This requires that
	 * a) the node has had transitive reduction performed on it (outgoing edges)
	 * and
	 * b) all parents have had transitive reduction performed (incoming edges)
	 * 
	 */
	private boolean isBlockingChildEmission(SgNode node) {
		switch (node.state) {
		case Removed:
		case UnderConstruction:
		case AwaitingReduction:
			assert(node.state != State.Removed);
			return true;
		case Reduced:
		case EmittedFromReduction:
		default:
			return false;
		}
	}
	private void load(SgNode node) {
		// Tell all parents that they are no longer waiting on us
		node.out.sort(SgEdge.BySequenceLength);
		node.state = SgNode.State.AwaitingReduction;
		node.mark = SgNode.TransitiveReductionMark.Vacant;
		node.waitingOnReduction = (int)node.out.stream().filter(n -> isBlockingParentTransitiveReduction(n.to)).count();
		node.waitingOnEmission = (int)node.in.stream().filter(n -> isBlockingChildEmission(n.from)).count();
		if (compressUnbranchingPaths && node.canCompress()) {
			SgEdge edge = new SgEdge(node);
			SgNode parent = edge.from;
			SgNode child = edge.to;
			if (!isBlockingChildEmission(parent)) {
				// we were blocking our child being emitted
				assert(parent.waitingOnEmission > 0);
				parent.waitingOnEmission--;
				if (parent.waitingOnEmission == 0 && parent.state == State.Reduced) {
					emit(parent);
				}
			}
			if (!isBlockingParentTransitiveReduction(child)) {
				// we were blocking TR because we we're loaded yet
				// the replacement child for our parent is now
				// not blocking TR
				assert(parent.waitingOnReduction > 0);
				parent.waitingOnReduction--;
				if (parent.waitingOnReduction == 0) {
					reduce(parent);
				}
			}
			return;
		}
		List<SgNode> toReduce = new ArrayList<>(node.in.size());
		for (SgEdge e : node.in) {
			SgNode parent = e.from;
			parent.waitingOnReduction--;
			if (parent.waitingOnReduction == 0) {
				toReduce.add(parent);
			}
		}
		for (SgNode parent : toReduce) {
			reduce(parent);
		}
		if (node.waitingOnReduction == 0) {
			reduce(node);
		}
	}
	/**
	 * Attempts to reduce the given node
	 * @param node
	 * @return
	 */
	private void reduce(SgNode node) {
		// perform transitive reduction
		int longest = 0;
		for (SgEdge vw : node.out) {
			SgNode w = vw.to;
			assert(vw.from == node);
			assert(w.state != State.UnderConstruction);
			w.mark = TransitiveReductionMark.InPlay;
			longest = Math.max(longest, vw.length());
		}
		longest += fuzz;
		
		for (SgEdge vw : node.out) {
			SgNode w = vw.to;
			if (w.mark == TransitiveReductionMark.InPlay) {
				for (SgEdge wx : w.out) {
					SgNode x = wx.to;
					if (wx.length() + vw.length() > longest) {
						// too long - no need to continue further with x
						break;
					}
					if (x.mark == TransitiveReductionMark.InPlay) {
						x.mark = TransitiveReductionMark.Eliminated;
					}
				}
			}
		}
		for (SgEdge vw : node.out) {
			SgNode w = vw.to;
			if (Defaults.SANITY_CHECK_ASSEMBLY_GRAPH && compressUnbranchingPaths) {
				// children we can performing transitive reduction over cannot be
				// compressable node or else the transitive reduction can fail
				assert(w.in.size() != 1 || w.out.size() != 1);
			}
			for (SgEdge wx : w.out) {
				SgNode x = wx.to;
				if (wx.length() < fuzz | wx.length() == w.out.get(0).length()) {
					if (x.mark == TransitiveReductionMark.InPlay) {
						x.mark = TransitiveReductionMark.Eliminated;
					}
				}
			}
		}
		List<SgNode> children = new ArrayList<>(node.out.size());
		List<SgEdge> toReduce = new ArrayList<>(node.out.size());
		for (SgEdge vw : node.out) {
			SgNode w = vw.to;
			children.add(w);
			if (w.mark == TransitiveReductionMark.Eliminated) {
				toReduce.add(vw);
			}
			w.mark = TransitiveReductionMark.Vacant;
		}
		// can't reduce the in above for loop since we're iterating over node.out
		for (SgEdge e : toReduce) {
			e.reduce();
		}
		node.state = State.Reduced;
		if (node.canCompress()) {
			// replace the node with the given edge
			SgEdge edge = new SgEdge(node);
			assert(edge.to.state != State.UnderConstruction);
			if (edge.to.state == State.AwaitingReduction) {
				// TODO: compression issues
				
			}
		} else {
			for (SgNode child : children) {
				child.waitingOnEmission--;
				// Tell our children (including those we removed edges to) that we have been reduced
				if (child.waitingOnEmission == 0 && child.state == State.Reduced) {
					emit(child);
				}
			}
			if (node.waitingOnEmission == 0) {
				output.add(node);
				node.state = State.EmittedFromReduction;
			}
		}
	}
	private void emit(SgNode node) {
		assert(node.state == State.Reduced);
		if (Defaults.SANITY_CHECK_ASSEMBLY_GRAPH) {
			assert(node.in.stream().allMatch(n -> !isBlockingChildEmission(node)));
		}
		output.add(node);
		node.state = State.EmittedFromReduction;
	}
	@Override
	public boolean hasNext() {
		ensureNext();
		return !output.isEmpty();
	}
	@Override
	public SgNode next() {
		ensureNext();
		if (output.isEmpty()) {
			throw new NoSuchElementException();
		}
		return output.poll();
	}
}