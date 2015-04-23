package au.edu.wehi.idsv.debruijn.positional;

import java.util.List;
import java.util.NavigableSet;
import java.util.Set;

public class Overview {
	public class Kmer { }
	public abstract class Evidence {
		Kmer[] kmers;
		int[] weight;
		String id;
		abstract int start();
		// first position is: [start - error, start + error]
		abstract int errorWidth(); 
		abstract int isAnchored(int offset);
	}
	// Data structures:
	public class UnanchoredEvidence extends Evidence {
		int lowStartPosition;
		int highStartPosition;
		int anchorPosition;
	}
	public class AnchoredEvidence extends Evidence {
		int startPosition;
		(int, int) anchoredOffsets;
	}
	Map<Kmer, IntervalTree<position, (offset, Evidence)> kmerMap;
	addEvidence(Evidence e) {
		for (int i = 0; i < e.size(); i++) {
			addToKmerMap(kmerMap, e.kmers[i],
				e.start() - e.errorWidth(),
				e.start() + e.errorWidth(),
				isAnchored(i), e);
		}
	}
	public class Node {
		Kmer kmer;
		int positionInterval;
		int weight;
		boolean isReference();
		List<(offset, Evidence)> getEvidence();
		PathNode containedIn;
		
		Node equivalent; // path collapsed
		// baseDiff(equivalent.kmer, kmer) < threshold (threshold applies over entire path)
		// equivalent.position = node.position
	}
	public class PathNode {
		List<Node> path; 
		// List<List<Node>> alternatePaths; // stored as node.equivalent
		// invariants
		// path[i+1].position = path[i].position + 1
		// alt[i][*].position = path[i].position
	}
	NavigableSet<PathNode> startingNodesByStartPosition;
	// Operations:
	// adds the given evidence to the graph
	add(evidence);
	// adds the given evidence kmer to the graph
	add(kmer, positionInterval, (weight, isReference, evidence));
	// performs bubble popping, leaf & path collapse on all PathNodes
	// entirely before the given position.
	simplify(before);
	// calls the best contig containing no nodes no or after the given position
	List<PathNode> callContigs(before)
	// remove all evidence contributing to the given contig from the graph
	removeEvidence(List<PathNode> contig) {
		HashSet<Evidence> support;
		// delete everything on path
		// delete any remaining support
	}
	removeEvidenceBefore(before);
}
