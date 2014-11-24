package au.edu.wehi.idsv.debruijn.subgraph;

import htsjdk.samtools.util.Log;

import java.io.File;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import au.edu.wehi.idsv.AssemblyEvidence;
import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.AssemblyFactory;
import au.edu.wehi.idsv.AssemblyParameters;
import au.edu.wehi.idsv.BreakendDirection;
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
			BreakendDirection direction,
			AssemblyParameters parameters,
			SubgraphAssemblyAlgorithmTrackerBEDWriter trackingWriter) {
		super(processContext, source, parameters.k, direction);
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
			return new SubgraphAlgorithmMetrics(processContext, referenceIndex, direction, ((TimedSubgraphSummary)ss).getCreationTime());
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
	private int subgraphMessageStartingSize = 4096;
	private boolean exceedsTimeout(SubgraphSummary ss) {
		return ss.isAnchored() && ss.getMaxAnchor() - ss.getMinAnchor() > getSafetyWidth();
	}
	private int getSafetyWidth() {
		return (int)(parameters.maxSubgraphFragmentWidth * source.getMaxConcordantFragmentSize());
	}
	private void visualisePrecollapsePathGraph(SubgraphSummary ss, PathGraphAssembler pga) {
		File directory = new File(new File(
				parameters.debruijnGraphVisualisationDirectory,
				processContext.getDictionary().getSequence(referenceIndex).getSequenceName()),
				String.format("%d", ss.getMinAnchor() - ss.getMinAnchor() % 100000));
		String filename = String.format("%s_%s_%d-%d_%d.precollapse.gexf",
				direction == BreakendDirection.Forward ? "f" : "b",
				processContext.getDictionary().getSequence(referenceIndex).getSequenceName(),
				ss.getMinAnchor(),
				ss.getMaxAnchor(),
				graphsExported); 
		directory.mkdirs();
		new StaticDeBruijnSubgraphPathGraphGexfExporter(this.parameters.k)
			.snapshot(pga)
			.saveTo(new File(directory, filename));
	}
	private void visualisePathGraph(SubgraphSummary ss, StaticDeBruijnSubgraphPathGraphGexfExporter graphExporter) {
		if (graphExporter == null) return;
		File directory = new File(new File(
				parameters.debruijnGraphVisualisationDirectory,
				processContext.getDictionary().getSequence(referenceIndex).getSequenceName()),
				String.format("%d", ss.getMinAnchor() - ss.getMinAnchor() % 100000));
		String filename = String.format("%s_%s_%d-%d_%d.subgraph.gexf",
				direction == BreakendDirection.Forward ? "f" : "b",
				processContext.getDictionary().getSequence(referenceIndex).getSequenceName(),
				ss.getMinAnchor(),
				ss.getMaxAnchor(),
				graphsExported++);
		directory.mkdirs();
		graphExporter.saveTo(new File(directory, filename));
	}
	private boolean shouldVisualise(boolean timeout) {
		return parameters.debruijnGraphVisualisationDirectory != null && (parameters.visualiseAll || (timeout && parameters.visualiseTimeouts));
	}
	/**
	 * Assembles contigs which do not have any relevance at or after the given position 
	 * @param position
	 * @return
	 */
	public Iterable<AssemblyEvidence> assembleContigsBefore(int position) {
		List<AssemblyEvidence> contigs = Lists.newArrayList();
		for (SubgraphSummary ss : subgraphs) {
			boolean timeoutExceeded = exceedsTimeout(ss);
			if (timeoutExceeded) {
				log.warn(String.format("Subgraph at %s:%d-%d exceeded maximum width of %d - calling",
						processContext.getDictionary().getSequence(this.referenceIndex).getSequenceName(),
						ss.getMinAnchor(),
						ss.getMaxAnchor(),
						getSafetyWidth()));
			}
			if ((ss.getMaxAnchor() < position || timeoutExceeded) && ss.isAnchored()) {
				SubgraphAssemblyAlgorithmTracker tracker = createSubgraphAssemblyAlgorithmTracker(ss);
				tracker.finalAnchors(ss.getMinAnchor(), ss.getMaxAnchor());
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
				int width = ss.getMaxAnchor() - ss.getMinAnchor();
				if (width > subgraphMessageStartingSize) {
					subgraphMessageStartingSize = width;
					log.debug(String.format("Subgraph width=%s [%d, %d] has %d paths: %s", width, ss.getMinAnchor(), ss.getMaxAnchor(), pga.getPathCount(), debugOutputKmerSpread(ss)));
				}
				StaticDeBruijnSubgraphPathGraphGexfExporter graphExporter = null;
				if (shouldVisualise(timeoutExceeded)) {
					graphExporter = new StaticDeBruijnSubgraphPathGraphGexfExporter(this.parameters.k);
				}
				for (List<Long> contig : pga.assembleContigs(graphExporter)) {
					AssemblyEvidence variant = toAssemblyEvidence(contig, tracker);
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
	private String debugOutputKmerSpread(SubgraphSummary ss) {
		long maxKmer = 0;
		int maxSpread = -1;
		for (long kmer : reachableFrom(ss.getAnyKmer())) {
			DeBruijnSubgraphNode node = getKmer(kmer);
			int refWidth = 0;
			if (node.getMinReferencePosition() != null) {
				refWidth = node.getMaxReferencePosition() - node.getMinReferencePosition();
			}
			int mateWidth = 0;
			if (node.getMinMatePosition() != null) {
				mateWidth = node.getMaxMatePosition() - node.getMinMatePosition();
			}
			if (Math.max(mateWidth, refWidth) > maxSpread) {
				maxSpread = Math.max(mateWidth, refWidth);
				maxKmer = kmer;
			}
		}
		return String.format("Max kmer %s ref:[%d,%d] mate:[%d,%d]",
				KmerEncodingHelper.toString(getK(), maxKmer),
				getKmer(maxKmer).getMinReferencePosition(),
				getKmer(maxKmer).getMaxReferencePosition(),
				getKmer(maxKmer).getMinMatePosition(),
				getKmer(maxKmer).getMaxMatePosition());
	}
	/**
	 * Removes all kmers not relevant at or after the given position
	 * @param position
	 */
	public void removeBefore(int position) {
		List<SubgraphSummary> toRemove = Lists.newArrayList();
		for (SubgraphSummary ss : subgraphs) {
			if (ss.getMaxAnchor() < position || exceedsTimeout(ss) || !ss.isAnchored()) {
				for (long kmer : reachableFrom(ss.getAnyKmer())) {
					remove(kmer);
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
	private AssemblyEvidence toAssemblyEvidence(List<Long> contigKmers, SubgraphAssemblyAlgorithmTracker tracker) {
		int refCount = -1;
		int refAnchor = 0;
		Integer mateAnchor = null;
		// Advance to first non-reference kmer
		while (getKmer(contigKmers.get(++refCount)).isReference());
		if (refCount > 0) {
			refAnchor = getKmer(contigKmers.get(refCount - 1)).getBestReferencePosition();
		}
		// Iterate over breakpoint kmers
		boolean messagePrinted = false;
		for (int i = refCount; i < contigKmers.size(); i++) {
			DeBruijnSubgraphNode kmer = getKmer(contigKmers.get(i));
			if (kmer.isReference()) continue;
			Integer mp = direction == BreakendDirection.Forward ? kmer.getMaxMatePosition() : kmer.getMinMatePosition();
			if (mateAnchor == null) {
				mateAnchor = mp;
			} else {
				if (mp != null) {
					// take closest mate anchor
					mateAnchor = direction == BreakendDirection.Forward ?
							Math.max(mp, mateAnchor) :
							Math.min(mp, mateAnchor);
				}
				if (mp == null && refCount == 0) {
					//         D-E    <-- we are attempting to assemble DE from soft-clips after assembling ABC 
					//        /
					//     A-B-C             (this shouldn't happen.)
					//    /
					// R-R
					// 
					//throw new RuntimeException("Sanity check failure: attempted to assemble anchored read as unanchored");
					
					// or we have multiple events at the same loci:
					//  RP - SC
					// =====>--------<===C
					//                 ==CCCCCC
					// read pair and long soft clip sharing same kmer
					// soft clip is sufficiently long that no anchor kmer exists
					if (!messagePrinted) {
						log.debug("Unsufficiently anchored SC being treated as mate evidence at ", kmer.getSupportingEvidenceList().iterator().next().getBreakendSummary().toString(processContext));
						messagePrinted = true;
					}
				}
			}
		}
		ContigAssembly ca = debruijnContigAssembly(contigKmers, refCount);
		AssemblyEvidence evidence = null;
		if (refCount > 0) {
			// anchored read
			evidence = AssemblyFactory.createAnchored(processContext, source, direction, ca.support,
					referenceIndex, refAnchor, refCount + getK() - 1,
					ca.baseCalls, ca.baseQuals, ca.normalBaseCount, ca.tumourBaseCount);
		} else if (mateAnchor != null) {
			// inexact breakend
			evidence = AssemblyFactory.createUnanchored(processContext, source, ca.support,
					ca.baseCalls, ca.baseQuals, ca.normalBaseCount, ca.tumourBaseCount);
		} else {
			// Assembly is neither anchored by breakend nor anchored by mate pair
			// This occurs when the de bruijn graph has a non-reference kmer fork
			// the best path will be taken first which removes the anchor from the
			// forked path. We are left with an unanchored assembly that we can't
			// do anything with
			//                      B-B-B   <- assembly B is unanchored
			//                     /
			//                A-A-A-A-A-A-A-A
			//               /
			// Ref: A-A-A-A-A
			//
			return null;
		}
		return evidence;
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
			if (node.getMinMatePosition() != null) assert(node.getSubgraph().getMinAnchor() <= node.getMinMatePosition());
			if (node.getMaxMatePosition() != null) assert(node.getSubgraph().getMaxAnchor() >= node.getMaxMatePosition());
			if (node.getMinReferencePosition() != null) assert(node.getSubgraph().getMinAnchor() <= node.getMinReferencePosition());
			if (node.getMaxReferencePosition() != null) assert(node.getSubgraph().getMaxAnchor() >= node.getMaxReferencePosition());
		}		
		//int lastMax = Integer.MIN_VALUE;
		for (SubgraphSummary ss : subgraphs) {
			//assert(ss.getMaxAnchor() >= lastMax); // sort order
			//lastMax = ss.getMaxAnchor();
			assert(ss == ss.getRoot()); // only root subgraphs should be in list
			//assert(ss.getRoot().getMinAnchor() != Integer.MAX_VALUE); // Will happen until we anchor long SC & reads split by Ns
			//assert(ss.getRoot().getMaxAnchor() != Integer.MIN_VALUE); // Will happen until we anchor long SC & reads split by Ns
		}
		return true;
	}
	public boolean sanityCheckSubgraphs(int minExpected, int maxExpected) {
		sanityCheckSubgraphs();
		for (SubgraphSummary ss : subgraphs) {
			assert(ss.getMaxAnchor() >= minExpected);
			assert(ss.getMaxAnchor() <= maxExpected);
		}
		return true;
	}
	public String getStateSummaryMetrics() {
		String result = String.format("kmers=%d subgraphs=%d", size(), subgraphs.size());
		if (size() > 0) {
			int minAnchor = Integer.MAX_VALUE;
			int maxAnchor = Integer.MIN_VALUE;
			int maxWidth = 0;
			for (SubgraphSummary ss : subgraphs) {
				minAnchor = Math.min(minAnchor, ss.getMinAnchor());
				maxAnchor = Math.max(maxAnchor, ss.getMaxAnchor());
				maxWidth = Math.max(maxWidth, ss.getMaxAnchor() - ss.getMinAnchor());
			}
			result = result + String.format(" maxWidth=%d anchor=[%d,%d]", maxWidth, minAnchor, maxAnchor);
		}
		return result;
	}
}
