package au.edu.wehi.idsv.debruijn.subgraph;

import gnu.trove.procedure.TLongObjectProcedure;
import htsjdk.samtools.util.Log;

import java.io.File;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import au.edu.wehi.idsv.AssemblyEvidence;
import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.AssemblyParameters;
import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.debruijn.DeBruijnVariantGraph;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.debruijn.ReadKmer;
import au.edu.wehi.idsv.debruijn.VariantEvidence;
import au.edu.wehi.idsv.util.AlgorithmRuntimeSafetyLimitExceededException;
import au.edu.wehi.idsv.visualisation.NontrackingSubgraphTracker;
import au.edu.wehi.idsv.visualisation.StaticDeBruijnSubgraphPathGraphGexfExporter;
import au.edu.wehi.idsv.visualisation.SubgraphAlgorithmMetrics;
import au.edu.wehi.idsv.visualisation.SubgraphAssemblyAlgorithmTracker;
import au.edu.wehi.idsv.visualisation.SubgraphAssemblyAlgorithmTrackerBEDWriter;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

public class DeBruijnReadGraph extends DeBruijnVariantGraph<DeBruijnSubgraphNode> {
	private static final Log log = Log.getInstance(DeBruijnReadGraph.class);
	private SubgraphAssemblyAlgorithmTrackerBEDWriter trackingWriter;
	/**
	 * Connected subgraphs
	 */
	private final Set<SubgraphSummary> subgraphs = Sets.newHashSet();
	private final int referenceIndex;
	private final AssemblyParameters parameters;
	private int graphsExported = 0;
	
	/**
	 * 
	 * @param k
	 * @param direction
	 * @param parameters 
	 * @param gexf 
	 */
	public DeBruijnReadGraph(
			ProcessingContext processContext,
			AssemblyEvidenceSource source,
			int referenceIndex,
			AssemblyParameters parameters,
			SubgraphAssemblyAlgorithmTrackerBEDWriter trackingWriter) {
		super(processContext, source, parameters.k);
		this.referenceIndex = referenceIndex;
		this.parameters = parameters;
		this.trackingWriter = trackingWriter;
	}
	@Override
	protected DeBruijnSubgraphNode createNode(VariantEvidence evidence, int readKmerOffset, ReadKmer kmer) {
		return new DeBruijnSubgraphNode(evidence, readKmerOffset, kmer);
	}
	/**
	 * Sets the subgraph for this new node
	 */
	@Override
	protected void onKmerAdded(long kmer, DeBruijnSubgraphNode node) {
		super.onKmerAdded(kmer, node);
		SubgraphSummary g = null;
		// check if we are connected to an existing subgraph
		for (long adjKmer : KmerEncodingHelper.adjacentStates(getK(), kmer)) {
			if (adjKmer == kmer) continue; // ignore loops back to ourself
			DeBruijnSubgraphNode adjNode = getKmer(adjKmer);
			if (adjNode != null) {
				SubgraphSummary gadj = adjNode.getSubgraph();
				if (g == null) {
					g = gadj;
				} else if (g != gadj) {
					// this kmer merges two subgraphs
					if (g.add(gadj)) {
						subgraphs.remove(gadj);
					}
				}
			}
		}
		if (g == null) {
			// nothing next to us -> new subgraph
			g = createSubgraphSummary(kmer);
			subgraphs.add(g);
		}
		node.setSubgraph(g);
		g.addNode(node);
		if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
			assert(sanityCheckSubgraphs());
		}
	}
	private SubgraphSummary createSubgraphSummary(long kmer) {
		if (parameters.trackAlgorithmProgress) {
			return new TimedSubgraphSummary(kmer);
		} else {
			return new SubgraphSummary(kmer);
		}
	}
	private SubgraphAssemblyAlgorithmTracker createSubgraphAssemblyAlgorithmTracker(SubgraphSummary ss) {
		if (parameters.trackAlgorithmProgress) {
			return new SubgraphAlgorithmMetrics(processContext, referenceIndex, ((TimedSubgraphSummary)ss).getCreationTime());
		} else {
			return new NontrackingSubgraphTracker();
		}
	}
	/**
	 * Updates subgraph aggregates based on the new evidence
	 */
	@Override
	public DeBruijnSubgraphNode add(long kmer, DeBruijnSubgraphNode node) {
		DeBruijnSubgraphNode result = super.add(kmer, node);
		// update subgraph bounds
		if (result.getSubgraph() != null) {
			result.getSubgraph().addNode(node);
		}
		return result;
	}
	private boolean exceedsTimeout(SubgraphSummary ss) {
		return ss.getMaxLinearPosition() - ss.getMinLinearPosition() > source.getAssemblyMaximumEvidenceDelay();
	}
	private void visualisePrecollapsePathGraph(SubgraphSummary ss, PathGraphAssembler pga) {
		File directory = new File(
				parameters.debruijnGraphVisualisationDirectory,
				processContext.getDictionary().getSequence(referenceIndex).getSequenceName());
		String filename = String.format("%d_%s.precollapse.gexf",
				graphsExported,
				processContext.getLinear().encodedIntervalToString(ss.getMinLinearPosition(), ss.getMaxLinearPosition()).replace("-", "_").replace(":", "_")); 
		directory.mkdirs();
		new StaticDeBruijnSubgraphPathGraphGexfExporter(this.parameters.k)
			.snapshot(pga)
			.saveTo(new File(directory, filename));
	}
	private void visualisePathGraph(SubgraphSummary ss, StaticDeBruijnSubgraphPathGraphGexfExporter graphExporter) {
		if (graphExporter == null) return;
		File directory = new File(
				parameters.debruijnGraphVisualisationDirectory,
				processContext.getDictionary().getSequence(referenceIndex).getSequenceName());
		String filename = String.format("%d_%s.subgraph.gexf",
				graphsExported++,
				processContext.getLinear().encodedIntervalToString(ss.getMinLinearPosition(), ss.getMaxLinearPosition()).replace("-", "_").replace(":", "_"));
		directory.mkdirs();
		graphExporter.saveTo(new File(directory, filename));
	}
	private boolean shouldVisualise(boolean timeout) {
		return parameters.debruijnGraphVisualisationDirectory != null && (parameters.visualiseAll || (timeout && parameters.visualiseTimeouts));
	}
	/**
	 * Assembles contigs which do not have any relevance at or after the given position 
	 * @param position linear position
	 * @return
	 */
	public Iterable<AssemblyEvidence> assembleContigsBefore(long position) {
		List<AssemblyEvidence> contigs = Lists.newArrayList();
		for (SubgraphSummary ss : subgraphs) {
			boolean timeoutExceeded = exceedsTimeout(ss);
			if (timeoutExceeded) {
				log.warn(String.format("Subgraph at %s:%d-%d exceeded maximum width of %d - calling",
						processContext.getDictionary().getSequence(this.referenceIndex).getSequenceName(),
						ss.getMinLinearPosition(),
						ss.getMaxLinearPosition(),
						source.getAssemblyMaximumEvidenceDelay()));
			}
			if (ss.getMaxLinearPosition() < position || timeoutExceeded) {
				SubgraphAssemblyAlgorithmTracker tracker = createSubgraphAssemblyAlgorithmTracker(ss);
				tracker.finalAnchors(processContext.getLinear().getReferenceIndex(ss.getMinLinearPosition()), processContext.getLinear().getReferenceIndex(ss.getMaxLinearPosition()));
				tracker.assemblyStarted();
				PathGraphAssembler pga = new PathGraphAssembler(this, this.parameters, ss.getAnyKmer(), tracker);
				if (parameters.maxBaseMismatchForCollapse > 0) {
					if (shouldVisualise(timeoutExceeded)) {
						visualisePrecollapsePathGraph(ss, pga);
					}
					try {
						// simplify graph
						pga.collapse(parameters.maxBaseMismatchForCollapse, parameters.collapseBubblesOnly);
					} catch (AlgorithmRuntimeSafetyLimitExceededException e) {
						timeoutExceeded = true;
					}
				}
				StaticDeBruijnSubgraphPathGraphGexfExporter graphExporter = null;
				if (shouldVisualise(timeoutExceeded)) {
					graphExporter = new StaticDeBruijnSubgraphPathGraphGexfExporter(this.parameters.k);
				}
				for (List<Long> contig : pga.assembleContigs(graphExporter)) {
					AssemblyEvidence variant = createAssembly(contig);
					if (variant != null) {
						contigs.add(variant);
					}
				}
				tracker.assemblyComplete();
				if (trackingWriter != null) {
					trackingWriter.write(tracker);
				}
				if (shouldVisualise(timeoutExceeded)) {
					visualisePathGraph(ss, graphExporter);
				}
			}
		}
		Collections.sort(contigs, DirectedEvidence.ByStartEnd);
		return contigs;
	}
	/**
	 * Removes all kmers not relevant at or after the given position
	 * @param position
	 */
	public void removeBefore(long position) {
		List<SubgraphSummary> toRemove = Lists.newArrayList();
		for (final SubgraphSummary ss : subgraphs) {
			if (ss.getMaxLinearPosition() < position || exceedsTimeout(ss)) {
				if (ss.getKmerCount() > getBackingStore().capacity() * 0.0625) {
					// TODO: do performance testing to find optimal bulk delete thresholds
					getBackingStore().retainEntries(new TLongObjectProcedure<DeBruijnSubgraphNode>() {
						@Override
						public boolean execute(long kmer, DeBruijnSubgraphNode node) {
							return node.getSubgraph() != ss;
						}
					});
				} else {
					for (long kmer : reachableFrom(ss.getAnyKmer())) {
						remove(kmer);
					}
				}
				toRemove.add(ss);
			}
		}
		subgraphs.removeAll(toRemove);
	}
	/**
	 * Only non-reference reads are considered supporting the breakpoint.
	 */
	@Override
	public Set<DirectedEvidence> getSupportingEvidence(Iterable<Long> path) {
		Set<DirectedEvidence> reads = Sets.newHashSet();
		for (Long kmer : path) {
			DeBruijnSubgraphNode node = getKmer(kmer);
			if (!node.isReference()) {
				reads.addAll(node.getSupportingEvidenceList());
			}
		}
		return reads;
	}
	private static Log expensiveSanityCheckLog = Log.getInstance(DeBruijnReadGraph.class);
	public boolean sanityCheckSubgraphs() {
		if (expensiveSanityCheckLog != null) {
			expensiveSanityCheckLog.warn("Expensive sanity checking is being performed. Performance will be poor");
			expensiveSanityCheckLog = null;
		}
		for (long kmer : getAllKmers()) {
			DeBruijnSubgraphNode node = getKmer(kmer);
			assert(node.getSubgraph() != null);
			assert(subgraphs.contains(node.getSubgraph()));
			assert(node.getSubgraph().getMaxLinearPosition() >= node.getMaxLinearPosition());
			assert(node.getSubgraph().getMinLinearPosition() <= node.getMinLinearPosition());
		}		
		//int lastMax = Integer.MIN_VALUE;
		for (SubgraphSummary ss : subgraphs) {
			//assert(ss.getMaxAnchor() >= lastMax); // sort order
			//lastMax = ss.getMaxAnchor();
			assert(ss == ss.getRoot()); // only root subgraphs should be in list
			//assert(ss.getRoot().getMinAnchor() != Long.MAX_VALUE); // Will happen until we anchor long SC & reads split by Ns
			//assert(ss.getRoot().getMaxAnchor() != Long.MIN_VALUE); // Will happen until we anchor long SC & reads split by Ns
		}
		return true;
	}
	public boolean sanityCheckSubgraphs(long minExpected, long maxExpected) {
		sanityCheckSubgraphs();
		for (SubgraphSummary ss : subgraphs) {
			assert(ss.getMaxLinearPosition() >= minExpected);
			assert(ss.getMaxLinearPosition() <= maxExpected);
		}
		return true;
	}
	public String getStateSummaryMetrics() {
		String result = String.format("kmers=%d subgraphs=%d", size(), subgraphs.size());
		if (size() > 0) {
			long minAnchor = Long.MIN_VALUE;
			long maxAnchor = Long.MAX_VALUE;
			int maxWidth = 0;
			for (SubgraphSummary ss : subgraphs) {
				minAnchor = Math.min(minAnchor, ss.getMinLinearPosition());
				maxAnchor = Math.max(maxAnchor, ss.getMaxLinearPosition());
				maxWidth = Math.max(maxWidth, (int)(ss.getMaxLinearPosition() - ss.getMinLinearPosition()));
			}
			result = result + String.format(" %s maxWidth=%d ", processContext.getLinear().encodedIntervalToString(minAnchor, maxAnchor), maxWidth);
		}
		return result;
	}
}
