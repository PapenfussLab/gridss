package scambler;

import java.util.ArrayDeque;
import java.util.Deque;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Set;

import scambler.SgNode.State;

/**
 * Streaming implementation of Myer 2005 String Graph transitive reduction algorithm
 * that also performs node compression.
 * 
 * Nodes are output when we can guarantee that no edges in or out of the node will
 * be modified after emission.
 * 
 * Transitive reduction can remove an edge coming into a node when any parent is processed
 * 
 * Transitive reduction can remove an edge coming out a node when the node itself is
 * 
 * Therefore, to emit, both the node and all parents must have had transitive reduction
 * performed on them.
 * 
 * To perform transitive reduction, both the node, and all children must have all
 * outgoing edges loaded.
 * 
 * 
 * TODO: node compression
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
	public StringGraphTransitiveCompressor(Iterator<SgNode> it) {
		this.it = it;
	}
	private void ensureNext() {
		while (it.hasNext() & output.isEmpty()) {
			load(it.next());
		}
	}
	private void load(SgNode node) {
		// Tell all parents that they are no longer waiting on us
		node.state = SgNode.State.AwaitingReduction;
		node.mark = SgNode.TransitiveReductionMark.Vacant;
		node.out.sort(SgEdge.BySequenceLength);
		node.waitingOnReduction = (int)node.out.stream().filter(n -> n.to.state == State.UnderConstruction).count();
		node.waitingOnEmission = (int)node.in.stream().filter(n -> n.to.state != State.Reduced).count();
		for (SgEdge e : node.in) {
			SgNode parent = e.from;
			parent.waitingOnReduction--;
			if (parent.waitingOnReduction == 0) {
				reduce(parent);
			}
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
		
		List<SgNode> children;
		// Tell our children (including those we removed edges to) that we have been reduced
		for (SgNode child : children) {
			child.waitingOnEmission--;
			if (child.waitingOnEmission == 0) {
				output.add(child);
			}
		}
		
		// Now we're waiting for our parents to have transitive reduction performed
		if (node.waitingOnEmission == 0) {
			output.add(node);
		}
	}
	private void eliminate(SgEdge node) {
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