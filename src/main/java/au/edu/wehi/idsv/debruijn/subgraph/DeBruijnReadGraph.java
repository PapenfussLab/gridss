package au.edu.wehi.idsv.debruijn.subgraph;

import htsjdk.samtools.util.Log;

import java.io.File;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import au.edu.wehi.idsv.AssemblyBuilder;
import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.AssemblyParameters;
import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.VariantContextDirectedEvidence;
import au.edu.wehi.idsv.debruijn.DeBruijnVariantGraph;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.debruijn.ReadKmer;
import au.edu.wehi.idsv.debruijn.VariantEvidence;
import au.edu.wehi.idsv.metrics.RelevantMetrics;
import au.edu.wehi.idsv.util.AlgorithmRuntimeSafetyLimitExceededException;
import au.edu.wehi.idsv.visualisation.StaticDeBruijnSubgraphPathGraphGexfExporter;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

public class DeBruijnReadGraph extends DeBruijnVariantGraph<DeBruijnSubgraphNode> {
	private final Log log = Log.getInstance(DeBruijnReadGraph.class);
	public static final String ASSEMBLER_NAME = "debruijn-s";
	/**
	 * Connected subgraphs
	 */
	private final Set<SubgraphSummary> subgraphs = Sets.newHashSet();
	private final int referenceIndex;
	private final AssemblyParameters parameters;
	private final RelevantMetrics metrics;
	private int graphsExported = 0;
	/**
	 * 
	 * @param k
	 * @param direction
	 * @param parameters 
	 * @param gexf 
	 */
	public DeBruijnReadGraph(ProcessingContext processContext, AssemblyEvidenceSource source, int referenceIndex, BreakendDirection direction, AssemblyParameters parameters) {
		this(processContext, source, referenceIndex, direction, parameters, null);
	}
	/**
	 * 
	 * @param k
	 * @param direction
	 * @param parameters 
	 * @param gexf 
	 */
	public DeBruijnReadGraph(ProcessingContext processContext, AssemblyEvidenceSource source, int referenceIndex, BreakendDirection direction, AssemblyParameters parameters, RelevantMetrics metrics) {
		super(processContext, source, parameters.k, direction);
		this.referenceIndex = referenceIndex;
		this.parameters = parameters;
		this.metrics = metrics;
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
					SubgraphSummary newg = SubgraphSummary.merge(g, gadj);
					if (newg == g) {
						subgraphs.remove(gadj);
					} else if (newg == gadj) {
						subgraphs.remove(g);
					} else {
						subgraphs.remove(gadj);
						subgraphs.remove(g);
						subgraphs.add(newg);
					}
					g = newg;
				}
			}
		}
		if (g == null) {
			// nothing next to us -> new subgraph
			g = new SubgraphSummary(kmer);
			subgraphs.add(g);
		}
		node.setSubgraph(g);
	}
	/**
	 * Updates subgraph aggregates based on the new evidence
	 */
	@Override
	public DeBruijnSubgraphNode add(long kmer, DeBruijnSubgraphNode node) {
		DeBruijnSubgraphNode result = super.add(kmer, node);
		// update subgraph bounds
		result.getSubgraph().addNode(node);
		//sanityCheckSubgraphs();
		return result;
	}
	private int maxSubgraphSize = 4096;
	private boolean exceedsTimeout(SubgraphSummary ss) {
		if (metrics != null) {
			return ss.getMaxAnchor() - ss.getMinAnchor() > parameters.maxSubgraphFragmentWidth * metrics.getMaxFragmentSize();
		}
		return false;
	}
	/**
	 * Assembles contigs which do not have any relevance at or after the given position 
	 * @param position
	 * @return
	 */
	public Iterable<VariantContextDirectedEvidence> assembleContigsBefore(int position) {
		if (position % 1000000 == 0) {
			log.debug(String.format("At %s:%d de Bruijn graph size=%d, subgraphs=%d", processContext.getDictionary().getSequence(this.referenceIndex).getSequenceName(), position, size(), subgraphs.size()));
		}
		List<VariantContextDirectedEvidence> contigs = Lists.newArrayList();
		for (SubgraphSummary ss : subgraphs) {
			boolean timeoutExceeded = exceedsTimeout(ss);
			if (timeoutExceeded) {
				log.warn(String.format("Subgraph at %s:%d-%d exceeded maximum width - calling immediately",
						processContext.getDictionary().getSequence(this.referenceIndex).getSequenceName(),
						ss.getMinAnchor(),
						ss.getMaxAnchor()));
			}
			if (ss.getMaxAnchor() < position || timeoutExceeded) {
				PathGraphAssembler pga = new PathGraphAssembler(this, this.parameters, ss.getAnyKmer());
				if (parameters.maxBaseMismatchForCollapse > 0) {
					if (parameters.debruijnGraphVisualisationDirectory != null && (parameters.visualiseAll || (timeoutExceeded && parameters.visualiseTimeouts))) {
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
					// simplify graph
					try {
						pga.collapse(parameters.maxBaseMismatchForCollapse, parameters.collapseBubblesOnly);
					} catch (AlgorithmRuntimeSafetyLimitExceededException e) {
						timeoutExceeded = true;
					}
				}
				int width = ss.getMaxAnchor() - ss.getMinAnchor();
				if (width > maxSubgraphSize) {
					maxSubgraphSize = width;
					log.debug(String.format("Subgraph width=%s [%d, %d] has %d paths: %s", width, ss.getMinAnchor(), ss.getMaxAnchor(), pga.getPaths().size(), debugOutputKmerSpread(ss)));
				}
				StaticDeBruijnSubgraphPathGraphGexfExporter graphExporter = null;
				if (parameters.debruijnGraphVisualisationDirectory != null && (parameters.visualiseAll || (timeoutExceeded && parameters.visualiseTimeouts))) {
					graphExporter = new StaticDeBruijnSubgraphPathGraphGexfExporter(this.parameters.k);
				}
				for (List<Long> contig : pga.assembleContigs(graphExporter)) {
					VariantContextDirectedEvidence variant = toAssemblyEvidence(contig);
					if (variant != null) {
						contigs.add(variant);
					}
				}
				if (graphExporter != null) {
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
			}
		}
		Collections.sort(contigs, VariantContextDirectedEvidence.ByLocationStart);
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
			if (ss.getMaxAnchor() < position || exceedsTimeout(ss)) {
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
				reads.addAll(node.getSupportingEvidence());
			}
		}
		return reads;
	}
	private VariantContextDirectedEvidence toAssemblyEvidence(List<Long> contigKmers) {
		int refCount = 0;
		int refAnchor = 0;
		Integer mateAnchor = null;
		// Iterate over reference anchor
		for (long kmer : contigKmers) {
			if (!getKmer(kmer).isReference()) break;
			refCount++;
			refAnchor = getKmer(kmer).getBestReferencePosition();
		}
		// Iterate over breakpoint kmers
		boolean messagePrinted = false;
		for (long kmer : contigKmers) {
			if (getKmer(kmer).isReference()) continue;
			Integer mp = direction == BreakendDirection.Forward ? getKmer(kmer).getMaxMatePosition() : getKmer(kmer).getMinMatePosition();
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
						log.debug("Unsufficiently anchored SC being treated as mate evidence at ", getKmer(kmer).getSupportingEvidence().iterator().next().getBreakendSummary().toString(processContext));
						messagePrinted = true;
					}
				}
			}
		}
		AssemblyBuilder builder = debruijnContigAssembly(contigKmers)
				.assemblerName(ASSEMBLER_NAME);
		if (refCount > 0) {
			// anchored read
			builder.referenceAnchor(referenceIndex, refAnchor)
				.anchorLength(refCount + getK() - 1);
		} else if (mateAnchor != null) {
			// inexact breakend
			builder.mateAnchor(referenceIndex, mateAnchor);
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
		return builder.makeVariant();
	}
	private static Log expensiveSanityCheckLog = Log.getInstance(DeBruijnReadGraph.class);
	public boolean sanityCheckSubgraphs() {
		if (expensiveSanityCheckLog != null) {
			expensiveSanityCheckLog.warn("Expensive sanity checking is being performed. Performance will be poor");
			expensiveSanityCheckLog = null;
		}
		for (long kmer : getAllKmers()) {
			DeBruijnSubgraphNode node = getKmer(kmer);
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
}
