package scambler;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.stream.Collectors;

import au.edu.wehi.idsv.Defaults;
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
 * Node compression:
 * node: >= loaded
 * 
 * Emission: (no changes to edges)
 * node: >= loaded
 * parent: loaded <= emitted
 * child: loaded <= emitted
 * 
 * 
 * 
 * @author Daniel Cameron
 *
 */
public class StringGraphTransitiveCompressor implements Iterator<SgNode> {
	private final Iterator<SgNode> it;
	private final Deque<SgNode> output = new ArrayDeque<>();
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
	private boolean canEmit(SgNode node) {
		return node.edgeState >= SgNode.EDGES_TRANSISTIVE_REDUCED &&
				node.in.stream().allMatch(n -> n.from.edgeState >= SgNode.EDGES_TRANSISTIVE_REDUCED) &&
				node.out.stream().allMatch(n -> n.to.edgeState >= SgNode.EDGES_TRANSISTIVE_REDUCED);
	}
	private boolean canReduce(SgNode node) {
		return node.edgeState >= SgNode.EDGES_CONSTRUCTED && node.edgeState < SgNode.EDGES_FINALISED &&
				node.out.stream().allMatch(n -> n.to.edgeState >= SgNode.EDGES_CONSTRUCTED && n.to.edgeState < SgNode.EDGES_FINALISED);
	}
	private boolean canCompress(SgNode node) {
		return node.edgeState >= SgNode.EDGES_CONSTRUCTED && node.edgeState < SgNode.EDGES_FINALISED &&
				node.in.stream().allMatch(n -> n.from.edgeState < SgNode.EDGES_FINALISED) &&
				node.out.stream().allMatch(n -> n.to.edgeState < SgNode.EDGES_FINALISED);
	}
	private boolean shouldReduce(SgNode node) {
		return node.edgeState == SgNode.EDGES_CONSTRUCTED && canReduce(node);
	}
	private boolean shouldCompress(SgNode node) {
		return compressUnbranchingPaths && node.isCompressible() && canCompress(node) &&
				(node.edgeState >= SgNode.EDGES_CONSTRUCTED || node.edgeState < SgNode.EDGES_FINALISED);
	}
	private boolean shouldEmit(SgNode node) {
		return node.edgeState == SgNode.EDGES_TRANSISTIVE_REDUCED &&
				canEmit(node);
	}
	private void load(SgNode node) {
		node.out.sort(SgEdge.BySequenceLength);
		node.mark = SgNode.TransitiveReductionMark.Vacant;
		node.edgeState = SgNode.EDGES_CONSTRUCTED;
		if (shouldCompress(node)) {
			SgEdge edge = new SgEdge(node);
			SgNode parent = edge.from;
			SgNode child = edge.to;
			assert(!shouldCompress(parent));
			assert(!shouldCompress(child));
			if (shouldReduce(parent)) {
				reduce(parent);
			}
			if (shouldEmit(child))  {
				emit(child);
			}
			// Loading a node will not change the reduction status of a child
			// since reduction of that node only depends on the node itself
			// and its children
			return;
		}
		List<SgNode> toReduce = new ArrayList<>(node.in.size());
		for (SgEdge e : node.in) {
			SgNode parent = e.from;
			if (parent.edgeState == SgNode.EDGES_CONSTRUCTED && canReduce(parent)) {
				toReduce.add(parent);
			}
		}
		for (SgNode parent : toReduce) {
			reduce(parent);
		}
		if (shouldReduce(node)) {
			reduce(node);
		}
	}
	/**
	 * Attempts to reduce the given node
	 * @param node
	 * @return
	 */
	private void reduce(SgNode node) {
		assert(shouldReduce(node));
		assert(node.edgeState == SgNode.EDGES_CONSTRUCTED);
		// perform transitive reduction
		int longest = 0;
		for (SgEdge vw : node.out) {
			SgNode w = vw.to;
			assert(vw.from == node);
			assert(w.edgeState >= SgNode.EDGES_CONSTRUCTED);
			assert(w.edgeState < SgNode.EDGES_FINALISED);
			w.mark = TransitiveReductionMark.InPlay;
			longest = Math.max(longest, vw.length());
		}
		longest += fuzz;
		
		for (SgEdge vw : node.out) {
			SgNode w = vw.to;
			if (w.mark == TransitiveReductionMark.InPlay) {
				for (SgEdge wx : w.out) {
					SgNode x = wx.to;
					assert(x.edgeState < SgNode.EDGES_FINALISED);
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
				// compressible node or else the transitive reduction can fail
				assert(!w.isCompressible());
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
		node.edgeState = SgNode.EDGES_TRANSISTIVE_REDUCED;
		// Compress nodes to retain invariant that
		// all no node has in and out degrees of 1
		if (compressUnbranchingPaths) {
			for (SgNode child : children) {
				if (shouldCompress(child)) {
					new SgEdge(child);
				}
			}
		}
		if (compressUnbranchingPaths && shouldCompress(node)) {
			SgEdge edge = new SgEdge(node);
			if (shouldReduce(edge.from)) {
			}
			if (shouldEmit(edge.to)) {
				emit(edge.to);
			}
		} else {
			// all parents might now be able to be reduced
			List<SgNode> parents = node.in.stream().map(e -> e.from).collect(Collectors.toList());
			for (SgNode parent : parents) {
				if (shouldReduce(parent)) {
					reduce(parent);
				}
				if (shouldEmit(parent)) {
					emit(parent);
				}
			}
			for (SgNode child : children) {
				if (shouldEmit(child)) {
					reduce(child);
				}
			}
		}
		if (shouldEmit(node)) {
			emit(node);
		}
	}
	private void emit(SgNode node) {
		assert(canEmit(node));
		assert(node.edgeState == SgNode.EDGES_TRANSISTIVE_REDUCED);
		node.edgeState = SgNode.EDGES_FINALISED;
		output.add(node);
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