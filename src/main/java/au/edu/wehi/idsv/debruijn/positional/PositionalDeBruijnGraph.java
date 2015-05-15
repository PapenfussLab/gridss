package au.edu.wehi.idsv.debruijn.positional;

import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.NavigableSet;
import java.util.TreeSet;

public class PositionalDeBruijnGraph {
	private Map<Long, NavigableSet<KmerSupportNode>> rawKmers = new HashMap<Long, NavigableSet<KmerSupportNode>>();
	public void addEvidence(Evidence e) {
		for (int i = 0; i < e.length(); i++) {
			add(e.node(i));
		}
	}
	private void add(KmerSupportNode node) {
		NavigableSet<KmerSupportNode> set = rawKmers.get(node.kmer());
		if (set == null) {
			set = new TreeSet<KmerSupportNode>();
			rawKmers.put(node.kmer(), set);
		}
		set.add(node);
	}
	public Map<Long, NavigableSet<KmerAggregateNode>> toKmerAggregateGraph() {
		Map<Long, NavigableSet<KmerAggregateNode>> result = new HashMap<Long, NavigableSet<KmerAggregateNode>>();
		for (Entry<Long, NavigableSet<KmerSupportNode>> support : rawKmers.entrySet()) {
			result.put(support.getKey(), KmerAggregateNode.aggregate(support.getKey(), support.getValue()));
		}
		return result;
	}
}
