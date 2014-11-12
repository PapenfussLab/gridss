package au.edu.wehi.idsv.debruijn.subgraph;

public class TimedSubgraphSummary extends SubgraphSummary { 
	private long startTime;
	public TimedSubgraphSummary(long startingKmer) {
		super(startingKmer);
		this.startTime = System.nanoTime();
	}
	@Override
	public boolean add(SubgraphSummary graph) {
		SubgraphSummary myRoot = graph.getRoot();
		SubgraphSummary graphRoot = graph.getRoot();
		if (myRoot instanceof TimedSubgraphSummary && graphRoot instanceof TimedSubgraphSummary) {
			((TimedSubgraphSummary)myRoot).startTime = Math.min(((TimedSubgraphSummary) myRoot).startTime, ((TimedSubgraphSummary) graphRoot).startTime);
		}
		return super.add(graph);
	}
	public long getCreationTime() {
		return startTime;
	}
}
