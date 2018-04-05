package au.edu.wehi.idsv.debruijn.positional;

import java.io.File;
import java.io.IOException;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.IdentityHashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Queue;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import com.google.common.collect.PeekingIterator;
import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;

import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.AssemblyFactory;
import au.edu.wehi.idsv.AssemblyIdGenerator;
import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.BreakendSummary;
import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.SanityCheckFailureException;
import au.edu.wehi.idsv.debruijn.DeBruijnGraphBase;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.graph.ScalingHelper;
import au.edu.wehi.idsv.model.Models;
import au.edu.wehi.idsv.sam.SamTags;
import au.edu.wehi.idsv.util.IntervalUtil;
import au.edu.wehi.idsv.util.MessageThrottler;
import au.edu.wehi.idsv.visualisation.AssemblyTelemetry.AssemblyChunkTelemetry;
import au.edu.wehi.idsv.visualisation.PositionalDeBruijnGraphTracker;
import au.edu.wehi.idsv.visualisation.PositionalDeBruijnGraphTracker.ContigStats;
import au.edu.wehi.idsv.visualisation.PositionalExporter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Log;
import it.unimi.dsi.fastutil.longs.Long2ObjectMap;
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.longs.LongOpenHashSet;
import it.unimi.dsi.fastutil.longs.LongSet;
import it.unimi.dsi.fastutil.objects.ObjectOpenCustomHashSet;


/**
 * Calls optimal contigs from a positional de Bruijn graph
 * 
 * @author Daniel Cameron
 *
 */
public class NonReferenceContigAssembler implements Iterator<SAMRecord> {
	private static final Log log = Log.getInstance(NonReferenceContigAssembler.class);
	/**
	 * Debugging tracker to ensure memoization export files have unique names
	 */
	private static final AtomicInteger pathExportCount = new AtomicInteger();
	private long telemetryLastflushContigs = System.nanoTime();
	private long telemetryLastflushReferenceNodes = System.nanoTime();
	private long telemetryLastloadGraphs = System.nanoTime();
	/**
	 * Since reference kmers are not scored, calculating 
	 * highest weighted results in a preference for paths
	 * ending at a RP with sequencing errors over a path
	 * anchored to the reference. 
	 * 
	 * To ensure that the anchored paths are scored higher
	 * than the unanchored paths, paths anchored to the
	 * reference are given a score adjustment larger than
	 * the largest expected score.
	 */
	static final int ANCHORED_SCORE = Integer.MAX_VALUE >> 2;
	/**
	 * TODO: check to see if this is worth doing
	 * Simplication reduces the graph size, but may trigger
	 * additional rememoization so might turn out to be a more
	 * expensive approach overall
	 */
	private static final boolean SIMPLIFY_AFTER_REMOVAL = false;
	private Long2ObjectMap<Collection<KmerPathNodeKmerNode>> graphByKmerNode = new Long2ObjectOpenHashMap<Collection<KmerPathNodeKmerNode>>();
	private TreeSet<KmerPathNode> graphByPosition = new TreeSet<KmerPathNode>(KmerNodeUtil.ByFirstStartKmer);
	private SortedSet<KmerPathNode> nonReferenceGraphByPosition = new TreeSet<KmerPathNode>(KmerNodeUtil.ByFirstStartKmer);
	private final EvidenceTracker evidenceTracker;
	private final AssemblyEvidenceSource aes;
	private final AssemblyIdGenerator assemblyNameGenerator;
	private final int maxEvidenceSupportIntervalWidth;
	private final int maxAnchorLength;
	private final int k;
	private final int referenceIndex;
	private final ContigStats stats = new ContigStats();
	private final PeekingIterator<KmerPathNode> underlying;
	private final String contigName;
	private final Queue<SAMRecord> called = new ArrayDeque<>();
	private int lastUnderlyingStartPosition = Integer.MIN_VALUE;
	private int lastNextPosition = Integer.MIN_VALUE;
	private RangeSet<Integer> toFlush = TreeRangeSet.create();
	private MemoizedContigCaller bestContigCaller;
	private int contigsCalled = 0;
	private long consumed = 0;
	private PositionalDeBruijnGraphTracker exportTracker = null;
	private AssemblyChunkTelemetry telemetry = null;
	public int getReferenceIndex() { return referenceIndex; }
	private int retainWidth() {
		return  maxContigAnchorLength() + Math.max(
				// safety check to ensure that flushed contigs don't call advanceUnderlying()
				maxExpectedBreakendLength() + minDistanceFromNextPositionForEvidenceToBeFullyLoaded() + maxAnchorLength + 1,
				// calculate retain width from contig 
				(int)(aes.getContext().getAssemblyParameters().positional.retainWidthMultiple * aes.getMaxConcordantFragmentSize())) - k + 1;
	}
	private int flushWidth() { return Math.max(1, (int)(aes.getContext().getAssemblyParameters().positional.flushWidthMultiple * aes.getMaxConcordantFragmentSize())) - k + 1; }
	private int maxExpectedBreakendLength() { return Math.max(2, ((int)(aes.getContext().getAssemblyParameters().maxExpectedBreakendLengthMultiple * aes.getMaxConcordantFragmentSize()) - k + 1)); }
	/**
	 * Worst case scenario is a RP providing single kmer support for contig
	 * read length - (k-1) + max-min fragment size
	 *
	 * ========== contig
	 *          --------- read contributes single kmer to contig  
	 *           \       \  in the earliest position
	 *            \  RP   \
	 *             \       \
	 *              \       \
	 *               ---------
	 *                       ^
	 *                       |
	 * Last position supported by this RP is here. 
	 */
	private int minDistanceFromNextPositionForEvidenceToBeFullyLoaded() {
		// TODO: work out why maxEvidenceSupportIntervalWidth + aes.getMaxReadLength() - k + 1 isn't sufficient distance
		return maxEvidenceSupportIntervalWidth + aes.getMaxReadLength() + aes.getMaxConcordantFragmentSize() + 1;
	}
	/**
	 * Creates a new structural variant positional de Bruijn graph contig assembly for the given chromosome
	 * @param it reads
	 * @param referenceIndex evidence source
	 * @param maxEvidenceDistance maximum distance from the first position of the first kmer of a read,
	 *  to the last position of the last kmer of a read. This should be set to read length plus
	 *  the max-min concordant fragment size
	 * @param maxAnchorLength maximum number of reference-supporting anchor bases to assemble
	 * @param k
	 * @param source assembly source
	 * @param assemblyNameGenerator 
	 * @param tracker evidence lookup
	 */
	public NonReferenceContigAssembler(
			Iterator<KmerPathNode> it,
			int referenceIndex,
			int maxEvidenceSupportIntervalWidth,
			int maxAnchorLength,
			int k,
			AssemblyEvidenceSource source,
			AssemblyIdGenerator assemblyNameGenerator,
			EvidenceTracker tracker,
			String contigName) {
		this.underlying = Iterators.peekingIterator(it);
		this.maxEvidenceSupportIntervalWidth = maxEvidenceSupportIntervalWidth;
		this.maxAnchorLength = maxAnchorLength;
		this.k = k;
		this.referenceIndex = referenceIndex;
		this.aes = source;
		this.assemblyNameGenerator = assemblyNameGenerator;
		this.evidenceTracker = tracker;
		this.contigName = contigName;
		initialiseBestCaller();
	}
	private void initialiseBestCaller() {
		this.bestContigCaller = new MemoizedContigCaller(ANCHORED_SCORE, maxEvidenceSupportIntervalWidth);
		for (KmerPathNode n : graphByPosition) {
			bestContigCaller.add(n);
		}
	}
	@Override
	public boolean hasNext() {
		ensureCalledContig();
		return !called.isEmpty();
	}
	@Override
	public SAMRecord next() {
		ensureCalledContig();
		return called.remove();
	}
	/**
	 * Flush calls outside retain window.
	 * This data
	 */
	private void flushCallsOutsideRetainWindow() {
		// Need to make sure our contig that we're forcing a call on is comprised of evidence that has
		// been fully loaded into the graph 
		//                    ------------------------------------------- contig
		//                    ^                                         ^                             
		// |<---flushWidth--->|<---------------------------------retainWidth------------------------------------>|
		// |             flushPosition                                                                     frontierStart
		// loadedStart                                                                                      nextPosition
		if (!nonReferenceGraphByPosition.isEmpty()) {
			int frontierStart = bestContigCaller.frontierStart(nextPosition());
			int flushPosition = frontierStart - retainWidth() - 1;
			int loadedStart = nonReferenceGraphByPosition.first().firstStart();
			if (loadedStart + flushWidth() < flushPosition) { // don't start flushing until we're at least flushWidth distance from the retain position
				ArrayDeque<KmerPathSubnode> forcedContig = null;
				// keep calling until we have no more contigs left even if we could be calling a suboptimal contig
				do {
					forcedContig = bestContigCaller.callBestContigStartingBefore(nextPosition(), flushPosition);
					if (forcedContig != null && forcedContig.getLast().lastEnd() + maxEvidenceSupportIntervalWidth >= nextPosition()) {
						// could happen if ap.removeMisassembledPartialContigsDuringAssembly is false as contig length would then be unbounded
						log.debug(String.format("Attempting to flush contig %s:%d-%d (+%d) when graph only loaded to %s:%d. Ignoring flush request.",
								contigName, forcedContig.getFirst().firstStart(), forcedContig.getLast().lastEnd(), maxEvidenceSupportIntervalWidth,
								contigName, nextPosition()));
						forcedContig = null;
					}
					callContig(forcedContig);
				} while (forcedContig != null);
				if (!called.isEmpty()) {
					//log.debug(String.format("Forced %d contigs in interval %s:%d-%d(%d)", called.size(), contigName, loadedStart, frontierStart, nextPosition()));
					if (getTelemetry() != null) {
						long currentTime = System.nanoTime();
						getTelemetry().flushContigs(referenceIndex, loadedStart, frontierStart, called.size(), currentTime - telemetryLastflushContigs);
						telemetryLastflushContigs = currentTime;
					}
					return;
				}
			}
		}
	}
	public int maxContigAnchorLength() {
		return maxContigAnchorLength(maxExpectedBreakendLength());
	}
	public int maxContigAnchorLength(int breakendContigLength) {
		return Math.max(breakendContigLength, maxAnchorLength);
	}
	private void ensureCalledContig() {
		while (called.isEmpty()) {
			flushExcessivelyDenseIntervals();
			// remove misassembled partial contigs
			if (aes.getContext().getAssemblyParameters().removeMisassembledPartialContigsDuringAssembly) {
				removeMisassembledPartialContig();
			}
			flushCallsOutsideRetainWindow(); // Safety calling to ensure loaded graph size is bounded
			flushReferenceNodes(); // Need to get rid of any orphaned reads that got error corrected to become fully reference-supporting
			if (!called.isEmpty()) {
				return;
			}
			// Call the next contig
			ArrayDeque<KmerPathSubnode> bestContig = bestContigCaller.bestContig(nextPosition());
			callContig(bestContig);
			if (called.isEmpty() && bestContig == null) {
				if (underlying.hasNext()) {
					advanceUnderlying();
				} else {
					flushReferenceNodes();
					if (!graphByPosition.isEmpty()) {
						if (!MessageThrottler.Current.shouldSupress(log, "non-empty graph with no contigs")) {
							log.error("Sanity check failure: non-empty graph with no contigs called " + contigName);
						}
					}
					return;
				}
			}
		}
		if (Defaults.SANITY_CHECK_MEMOIZATION) {
			assert(bestContigCaller.sanityCheckFrontier(nextPosition()));
			verifyMemoization();
		}
	}
	private int flushReferenceNodes_debug_message_count = 0;
	private void flushReferenceNodes() {
		int endPosition = nonReferenceGraphByPosition.isEmpty() ? nextPosition() : nonReferenceGraphByPosition.first().firstStart();
		// first position at which we are guaranteed to not be involved in any contig anchor sequence
		endPosition -= minDistanceFromNextPositionForEvidenceToBeFullyLoaded() + maxContigAnchorLength();
		if (!graphByPosition.isEmpty() && graphByPosition.first().lastEnd() < endPosition) {
			int startPosition = graphByPosition.first().firstStart();
			Collection<KmerPathSubnode> nodes = new ArrayList<>();
			for (KmerPathNode tn : graphByPosition) {
				if (tn.lastEnd() >= endPosition) {
					break;
				}
				if (tn.isReference()) {
					nodes.add(new KmerPathSubnode(tn));
				} else {
					if (flushReferenceNodes_debug_message_count == 0) {
						log.debug(String.format("Sanity check failure when flushing reference nodes before %s:%d. Found non-reference node starting at %d", contigName, endPosition, tn.firstStart()));
						flushReferenceNodes_debug_message_count++;
						break;
					}
				}
			}
			Set<KmerEvidence> toRemove = evidenceTracker.untrack(nodes);
			removeFromGraph(toRemove);
			if (getTelemetry() != null) {
				long currentTime = System.nanoTime();
				getTelemetry().flushReferenceNodes(referenceIndex, startPosition, endPosition, toRemove.size(), currentTime - telemetryLastflushReferenceNodes);
				telemetryLastflushReferenceNodes = currentTime;
			}
		}
	}
	/**
	 * Removes partial contigs that are longer than the maximum theoretical breakend contig length
	 */
	private void removeMisassembledPartialContig() {
		// |<---      maxExpectedBreakendLength            --->|
		//            |<--- maxEvidenceSupportIntervalWidth--->|   
		//                                                  nextPosition
		ArrayDeque<KmerPathSubnode> misassembly = bestContigCaller.frontierPath(nextPosition(), nextPosition() - maxExpectedBreakendLength());
		if (misassembly == null) return;
		// To be sure that all reads on the contig to remove have
		// been fully loaded, we don't remove nodes that could contain
		// a read that also contributed to an unprocessed node
		List<KmerPathSubnode> misassemblyToRemove = misassembly.stream()
			.filter(sn -> sn.lastEnd() + minDistanceFromNextPositionForEvidenceToBeFullyLoaded() < nextPosition())
			.collect(Collectors.toList());
		if (misassemblyToRemove.size() == 0) {
			return;
		}
		Set<KmerEvidence> evidence = evidenceTracker.untrack(misassemblyToRemove);
		removeFromGraph(evidence);
	}
	private int nextPosition() {
		if (!underlying.hasNext()) return Integer.MAX_VALUE;
		return underlying.peek().firstStart();
	}
	/**
	 * Loads additional nodes into the graph
	 * 
	 * By loaded in batches, we reduce our memoization frontier advancement overhead
	 */
	private void advanceUnderlying() {
		int loadUntil = nextPosition();
		if (loadUntil < Integer.MAX_VALUE) {
			loadUntil += minDistanceFromNextPositionForEvidenceToBeFullyLoaded();
		}
		advanceUnderlying(loadUntil);
	}
	/**
	 * Ensures that portions of the graph exceeding the maximum density are flushed
	 */
	private void flushExcessivelyDenseIntervals() {
		while (!toFlush.isEmpty()) {
			if (graphByPosition.isEmpty()) break;
			Range<Integer> range = toFlush.asRanges().iterator().next();
			toFlush.remove(range);
			int flushOnOrAfter = range.lowerEndpoint();
			int flushBefore = range.upperEndpoint();
			int advanceTo = flushBefore + minDistanceFromNextPositionForEvidenceToBeFullyLoaded();
			advanceUnderlying(advanceTo);
			List<KmerPathSubnode> toRemove = new ArrayList<>();
			Iterator<KmerPathNode> it = graphByPosition.descendingIterator(); 
			while (it.hasNext()) {
				KmerPathNode pn = it.next();
				if (pn.firstStart() < flushOnOrAfter) break;
				if (pn.firstStart() >= flushBefore) continue;
				toRemove.add(new KmerPathSubnode(pn));
			}
			Set<KmerEvidence> evidenceToRemove = evidenceTracker.untrack(toRemove);
			if (!evidenceToRemove.isEmpty()) { // it could all overlap our previous flush range
				removeFromGraph(evidenceToRemove);
			}
		}
	}
	/**
	 * Advances the graph to the given position
	 * @param loadUntil
	 * @return true if the nodes loaded exceed the maximum density and should be flushed
	 */
	private void advanceUnderlying(int loadUntil) {
		if (loadUntil < nextPosition()) return;
		lastNextPosition = nextPosition();
		int count = 0;
		while (underlying.hasNext() && nextPosition() <= loadUntil) {
			KmerPathNode node = underlying.next();
			assert(lastUnderlyingStartPosition <= node.firstStart());
			lastUnderlyingStartPosition = node.firstStart();
			if (Defaults.SANITY_CHECK_ASSEMBLY_GRAPH) {
				assert(evidenceTracker.matchesExpected(new KmerPathSubnode(node)));
			}
			addToGraph(node);
			consumed++;
			if (!node.isReference()) {
				count++;
			}
		}
		int advanceWidth = nextPosition() - lastNextPosition;
		float density = advanceWidth <= 0 ? 0 : count / (float)advanceWidth;
		boolean filtered = false;
		if (density > aes.getContext().getAssemblyParameters().positional.maximumNodeDensity) {
			//log.debug(String.format("Density of %.2f at %s:%d-%d exceeds maximum: excluding from assembling.", density, contigName, lastNextPosition, nextPosition()));
			toFlush.add(Range.closedOpen(lastNextPosition, nextPosition()));
			filtered = true;
		}
		if (getTelemetry() != null) {
			long currentTime = System.nanoTime();
			getTelemetry().loadGraph(referenceIndex, lastNextPosition, nextPosition(), count, filtered, currentTime - telemetryLastloadGraphs);
			telemetryLastloadGraphs = currentTime;
		}
	}
	/**
	 * Verifies that the memoization matches a freshly calculated memoization 
	 * @param contig
	 */
	private boolean verifyMemoization() {
		int preGraphSize = graphByPosition.size();
		MemoizedContigCaller mcc = new MemoizedContigCaller(ANCHORED_SCORE, maxEvidenceSupportIntervalWidth);
		for (KmerPathNode n : graphByPosition) {
			mcc.add(n);
		}
		mcc.bestContig(nextPosition());
		bestContigCaller.sanityCheckMatches(mcc);
		int postGraphSize = graphByPosition.size();
		assert(preGraphSize == postGraphSize);
		return true;
	}
	/**
	 * Calls the contig.
	 * Note: to call the contig, additional graph nodes may be required to be loaded. No loading checks
	 * (such as flushing excessively dense intervals) are performed on these additional loaded nodes.
	 * @param contig to call 
	 * @return SAMRecord containing breakend assembly, null if a valid break-end assembly contig
	 * could not be created from this contig 
	 */
	private SAMRecord callContig(ArrayDeque<KmerPathSubnode> rawcontig) {
		if (rawcontig == null) return null;
		ArrayDeque<KmerPathSubnode> contig = rawcontig;
		if (containsKmerRepeat(contig)) {
			// recalculate the called contig, this may break the contig at the repeated kmer
			MisassemblyFixer fixed = new MisassemblyFixer(contig);
			contig = new ArrayDeque<KmerPathSubnode>(fixed.correctMisassignedEvidence(evidenceTracker.support(contig)));
		}
		if (contig.isEmpty()) return null;
		
		int contigLength = contig.stream().mapToInt(sn -> sn.length()).sum();
		int targetAnchorLength = Math.max(Math.min(contigLength, maxExpectedBreakendLength()), maxAnchorLength);
		KmerPathNodePath startAnchorPath = new KmerPathNodePath(contig.getFirst(), false, targetAnchorLength + maxEvidenceSupportIntervalWidth + contig.getFirst().length());
		startAnchorPath.greedyTraverse(true, false);
		ArrayDeque<KmerPathSubnode> startingAnchor = startAnchorPath.headNode().asSubnodes();
		startingAnchor.removeLast();
		// make sure we have enough of the graph loaded so that when
		// we traverse forward, our anchor sequence will be fully defined
		advanceUnderlying(contig.getLast().lastEnd() + targetAnchorLength + minDistanceFromNextPositionForEvidenceToBeFullyLoaded()); 
		
		KmerPathNodePath endAnchorPath = new KmerPathNodePath(contig.getLast(), true, targetAnchorLength + maxEvidenceSupportIntervalWidth + contig.getLast().length());
		endAnchorPath.greedyTraverse(true, false);
		ArrayDeque<KmerPathSubnode> endingAnchor = endAnchorPath.headNode().asSubnodes();
		endingAnchor.removeFirst();
		
		List<KmerPathSubnode> fullContig = new ArrayList<KmerPathSubnode>(contig.size() + startingAnchor.size() + endingAnchor.size());
		fullContig.addAll(startingAnchor);
		fullContig.addAll(contig);
		fullContig.addAll(endingAnchor);
		
		byte[] bases = KmerEncodingHelper.baseCalls(fullContig.stream().flatMap(sn -> sn.node().pathKmers().stream()).collect(Collectors.toList()), k);
		byte[] quals = DeBruijnGraphBase.kmerWeightsToBaseQuals(k, fullContig.stream().flatMapToInt(sn -> sn.node().pathWeights().stream().mapToInt(Integer::intValue)).toArray());
		assert(quals.length == bases.length);
		// left aligned anchor position although it shouldn't matter since anchoring should be a single base wide
		int startAnchorPosition = startingAnchor.size() == 0 ? 0 : startingAnchor.getLast().lastStart() + k - 1;
		int endAnchorPosition = endingAnchor.size() == 0 ? 0 : endingAnchor.getFirst().firstStart();
		int startAnchorBaseCount = startingAnchor.size() == 0 ? 0 : startingAnchor.stream().mapToInt(n -> n.length()).sum() + k - 1;
		int endAnchorBaseCount = endingAnchor.size() == 0 ? 0 : endingAnchor.stream().mapToInt(n -> n.length()).sum() + k - 1;
		int startBasesToTrim = Math.max(0, startAnchorBaseCount - targetAnchorLength);
		int endingBasesToTrim = Math.max(0, endAnchorBaseCount - targetAnchorLength);
		bases = Arrays.copyOfRange(bases, startBasesToTrim, bases.length - endingBasesToTrim);
		quals = Arrays.copyOfRange(quals, startBasesToTrim, quals.length - endingBasesToTrim);
		
		Set<KmerEvidence> evidence = evidenceTracker.untrack(contig);
		List<DirectedEvidence> evidenceIds = evidence.stream().map(e -> e.evidence()).collect(Collectors.toList());
		SAMRecord assembledContig;
		if (startingAnchor.size() == 0 && endingAnchor.size() == 0) {
			assert(startBasesToTrim == 0);
			assert(endingBasesToTrim == 0);
			// unanchored
			if (evidence.size() > 0) {
				BreakendSummary be = Models.calculateBreakend(aes.getContext().getLinear(),
					evidence.stream().map(e -> e.evidence().getBreakendSummary()).collect(Collectors.toList()),
					evidence.stream().map(e -> ScalingHelper.toScaledWeight(e.evidenceQuality())).collect(Collectors.toList()));
				assembledContig = AssemblyFactory.createUnanchoredBreakend(aes.getContext(), aes, assemblyNameGenerator,
						be,
						evidenceIds,
						bases, quals);
				if (evidence.stream().anyMatch(e -> e.isAnchored())) {
					log.debug(String.format("Unanchored assembly %s at %s:%d contains anchored evidence", assembledContig.getReadName(), contigName, contig.getFirst().firstStart()));
				}
			} else {
				// Can't call contig because we don't known the supporting breakend positions 
				assembledContig = null;
			}
		} else if (startingAnchor.size() == 0) {
			// end anchored
			assembledContig = AssemblyFactory.createAnchoredBreakend(aes.getContext(), aes, assemblyNameGenerator,
					BreakendDirection.Backward, evidenceIds,
					referenceIndex, endAnchorPosition, endAnchorBaseCount - endingBasesToTrim,
					bases, quals);
		} else if (endingAnchor.size() == 0) {
			// start anchored
			assembledContig = AssemblyFactory.createAnchoredBreakend(aes.getContext(), aes, assemblyNameGenerator,
					BreakendDirection.Forward, evidenceIds,
					referenceIndex, startAnchorPosition, startAnchorBaseCount - startBasesToTrim,
					bases, quals);
		} else {
			if (startAnchorBaseCount + endAnchorBaseCount >= quals.length) {
				// no unanchored bases - not an SV assembly
				assembledContig = null;
			} else {
				assembledContig = AssemblyFactory.createAnchoredBreakpoint(aes.getContext(), aes, assemblyNameGenerator, evidenceIds,
						referenceIndex, startAnchorPosition, startAnchorBaseCount - startBasesToTrim,
						referenceIndex, endAnchorPosition, endAnchorBaseCount - endingBasesToTrim,
						bases, quals);
			}
		}
		if (contigLength > maxExpectedBreakendLength()) {
			log.debug(String.format("Called breakend contig %s at %s:%d-%d of length %d when maximum expected breakend contig length is %d",
					assembledContig == null ? "" : assembledContig.getReadName(),
					contigName,
					rawcontig.getFirst().firstStart(),
					rawcontig.getLast().lastEnd(),
					contigLength,
					maxExpectedBreakendLength()));
		}
		if (assembledContig != null) {
			if (aes.getContext().getConfig().getVisualisation().assemblyGraph) {
				try {
					PositionalExporter.exportDot(new File(aes.getContext().getConfig().getVisualisation().directory, "assembly." + contigName + "." + assembledContig.getReadName() + ".dot"), k, graphByPosition, fullContig);
				} catch (Exception ex) {
					log.debug(ex, "Error exporting assembly ", assembledContig != null ? assembledContig.getReadName() : "(null)", " ", contigName);
				}
			}
			if (aes.getContext().getConfig().getVisualisation().assemblyGraphFullSize) {
				try {
					PositionalExporter.exportNodeDot(new File(aes.getContext().getConfig().getVisualisation().directory, "assembly.fullsize." + contigName + "." + assembledContig.getReadName() + ".dot"), k, graphByPosition, fullContig);
				} catch (Exception ex) {
					log.debug(ex, "Error exporting assembly ", assembledContig != null ? assembledContig.getReadName() : "(null)", " ", contigName);
				}
			}
			if (aes.getContext().getConfig().getVisualisation().assemblyContigMemoization) {
				File file = new File(aes.getContext().getConfig().getVisualisation().directory, "assembly.path.memoization." + contigName + "." + Integer.toString(pathExportCount.incrementAndGet()) + ".csv");
				try {
					bestContigCaller.exportState(file);
				} catch (IOException e) {
					log.debug(e, " Unable to export assembly path memoization to ", file, " ", contigName);
				}
			}
		}
		stats.contigNodes = contig.size();
		stats.truncatedNodes = rawcontig.size() - contig.size();
		stats.contigStartPosition = contig.getFirst().firstStart();
		stats.startAnchorNodes = startingAnchor.size();
		stats.endAnchorNodes = endingAnchor.size();
		if (exportTracker != null) {
			exportTracker.trackAssembly(bestContigCaller);
		}
		if (assembledContig != null) {
			String categorySupport = PositionalContigCategorySupportHelper.getCategorySupport(fullContig, evidence, k);
			assembledContig.setAttribute(SamTags.ASSEMBLY_CATEGORY_COVERAGE_CIGAR, categorySupport);
		}
		// remove all evidence contributing to this assembly from the graph
		if (evidence.size() > 0) {
			removeFromGraph(evidence);
			if (Defaults.SANITY_CHECK_MEMOIZATION) {
				bestContigCaller.sanityCheck(graphByPosition);
			}
		} else {
			if (!MessageThrottler.Current.shouldSupress(log, "unsupported paths")) {
				log.error("Sanity check failure: found path with no support. Attempting to recover by direct node removal ", contigName);
			}
			for (KmerPathSubnode n : contig) {
				removeFromGraph(n.node(), true);
			}
		}
		contigsCalled++;
		if (assembledContig != null) {
			called.add(assembledContig);
		}
		return assembledContig;
	}
	private boolean containsKmerRepeat(Collection<KmerPathSubnode> contig) {
		LongSet existing = new LongOpenHashSet();
		for (KmerPathSubnode n : contig) {
			for (int i = 0; i < n.length(); i++) {
				if (!existing.add(n.node().kmer(i))) {
					return true;
				}
			}
			for (long kmer : n.node().collapsedKmers()) {
				if (!existing.add(kmer)) {
					return true;
				}
			}
		}
		return false;
	}
	/**
	 * Removes all evidence from the current graph
	 * @param evidence
	 */
	private void removeFromGraph(Set<KmerEvidence> evidence) {
		assert(!evidence.isEmpty());
		// Tracks what we need to remove from each kmer of each path node
		Map<KmerPathNode, List<List<KmerNode>>> toRemove = new IdentityHashMap<KmerPathNode, List<List<KmerNode>>>();
		for (KmerEvidence e : evidence) {
			updateRemovalList(toRemove, e);
		}
		if (toRemove.size() > aes.getContext().getAssemblyParameters().positional.forceFullMemoizationRecalculationAt * graphByPosition.size()) {
			bestContigCaller = null;
		}
		if (bestContigCaller != null) {
			// removes all KmerPathNodes that need mutation from the memoization 
			bestContigCaller.remove(toRemove.keySet());
		}
		Set<KmerPathNode> simplifyCandidates = null;
		if (SIMPLIFY_AFTER_REMOVAL) {
			simplifyCandidates = new ObjectOpenCustomHashSet<KmerPathNode>(new KmerPathNode.HashByFirstKmerStartPositionKmer<KmerPathNode>());
		}
		for (Entry<KmerPathNode, List<List<KmerNode>>> entry : toRemove.entrySet()) {
			// removing down-weighted replacement nodes from memoization is unnecessary since we just
			// removed the entire node earlier in this function
			removeWeight(entry.getKey(), entry.getValue(), simplifyCandidates, false);
		}
		if (SIMPLIFY_AFTER_REMOVAL) {
			simplify(simplifyCandidates);
		}
		if (Defaults.SANITY_CHECK_ASSEMBLY_GRAPH) {
			assert(sanityCheck());
			assert(sanityCheckDisjointNodeIntervals());
		}
		if (Defaults.SANITY_CHECK_MEMOIZATION && bestContigCaller != null) {
			// Force memoization recalculation now
			bestContigCaller.bestContig(nextPosition());
			// so we can check that our removal was correct
			verifyMemoization();
		}
		if (bestContigCaller == null) {
			initialiseBestCaller();
		}
	}
	/**
	 * Attempts to simplify the given nodes
	 * @param simplifyCandidates
	 */
	private void simplify(Set<KmerPathNode> simplifyCandidates) {
		while (!simplifyCandidates.isEmpty()) {
			simplify(simplifyCandidates.iterator().next(), simplifyCandidates);
		}
	}
	private void simplify(KmerPathNode node, Set<KmerPathNode> simplifyCandidates) {
		simplifyCandidates.remove(node);
		if (node.lastEnd() >= nextPosition() - 1) {
			// don't simplify graph if we haven't actually loaded all the relevant nodes
			return;
		}
		KmerPathNode prev = node.prevToMergeWith();
		if (prev != null && prev.lastEnd() < nextPosition() - 1) {
			simplifyCandidates.remove(prev);
			removeFromGraph(node, true);
			removeFromGraph(prev, true);
			node.prepend(prev);
			addToGraph(node);
		}
		KmerPathNode next = node.nextToMergeWith();
		if (next != null && next.lastEnd() < nextPosition() - 1) {
			simplifyCandidates.remove(next);
			removeFromGraph(node, true);
			removeFromGraph(next, true);
			next.prepend(node);
			addToGraph(next);
		}
	}
	private void updateRemovalList(Map<KmerPathNode, List<List<KmerNode>>> toRemove, KmerEvidence e) {
		for (int i = e.length() - 1; i >= 0; i--) {
			KmerSupportNode support = e.node(i);
			if (support != null) {
				if (support.lastEnd() >= nextPosition()) {
					String msg = String.format("Sanity check failure: %s extending to %d removed when input at %s:%d", e, support.lastEnd(), contigName, nextPosition());
					log.error(msg);
					throw new SanityCheckFailureException(msg);
				}
				updateRemovalList(toRemove, support);
			}
		}
	}
	private void updateRemovalList(Map<KmerPathNode, List<List<KmerNode>>> toRemove, KmerSupportNode support) {
		Collection<KmerPathNodeKmerNode> kpnknList = graphByKmerNode.get(support.lastKmer());
		if (kpnknList != null) {
			for (KmerPathNodeKmerNode n : kpnknList) {
				if (IntervalUtil.overlapsClosed(support.lastStart(), support.lastEnd(), n.lastStart(), n.lastEnd())) {
					updateRemovalList(toRemove, n, support);
				}
			}
		}
	}
	private void updateRemovalList(Map<KmerPathNode, List<List<KmerNode>>> toRemove, KmerPathNodeKmerNode node, KmerSupportNode support) {
		KmerPathNode pn = node.node();
		List<List<KmerNode>> list = toRemove.get(pn);
		if (list == null) {
			list = new ArrayList<List<KmerNode>>(pn.length());
			toRemove.put(pn, list);
		}
		int offset = node.offsetOfPrimaryKmer();
		while (list.size() <= offset) {
			list.add(null);
		}
		List<KmerNode> evidenceList = list.get(offset); 
		if (evidenceList == null) {
			evidenceList = new ArrayList<KmerNode>();
			list.set(offset, evidenceList);
		}
		evidenceList.add(support);
	}
	private void removeWeight(KmerPathNode node, List<List<KmerNode>> toRemove, Set<KmerPathNode> simplifyCandidates, boolean includeMemoizationRemoval) {
		if (node == null) return;
		assert(node.length() >= toRemove.size());
		// remove from graph
		removeFromGraph(node, includeMemoizationRemoval);
		if (SIMPLIFY_AFTER_REMOVAL) {
			simplifyCandidates.addAll(node.next());
			simplifyCandidates.addAll(node.prev());
			simplifyCandidates.remove(node);
		}
		Collection<KmerPathNode> replacementNodes = KmerPathNode.removeWeight(node, toRemove);
		for (KmerPathNode split : replacementNodes) {
			if (Defaults.SANITY_CHECK_ASSEMBLY_GRAPH) {
				assert(evidenceTracker.matchesExpected(new KmerPathSubnode(split)));
			}
			addToGraph(split);
		}
		if (SIMPLIFY_AFTER_REMOVAL) {
			simplifyCandidates.addAll(replacementNodes);
		}
	}
	private void addToGraph(KmerPathNode node) {
		boolean added = graphByPosition.add(node);
		assert(added);
		if (!node.isReference()) {
			nonReferenceGraphByPosition.add(node);
		}
		for (int i = 0; i < node.length(); i++) {
			addToGraph(new KmerPathNodeKmerNode(node, i));
		}
		for (int i = 0; i < node.collapsedKmers().size(); i++) {
			addToGraph(new KmerPathNodeKmerNode(i, node));
		}
		if (bestContigCaller != null) {
			bestContigCaller.add(node);
		}
	}
	private void removeFromGraph(KmerPathNode node, boolean includeMemoizationRemoval) {
		if (includeMemoizationRemoval) {
			if (bestContigCaller != null) {
				bestContigCaller.remove(node);
			}
		}
		boolean removed = graphByPosition.remove(node);
		nonReferenceGraphByPosition.remove(node);
		assert(removed);
		for (int i = 0; i < node.length(); i++) {
			removeFromGraph(new KmerPathNodeKmerNode(node, i));
		}
		for (int i = 0; i < node.collapsedKmers().size(); i++) {
			removeFromGraph(new KmerPathNodeKmerNode(i, node));
		}
	}
	private void addToGraph(KmerPathNodeKmerNode node) {
		Collection<KmerPathNodeKmerNode> list = graphByKmerNode.get(node.firstKmer());
		if (list == null) {
			list = new ArrayList<KmerPathNodeKmerNode>();
			graphByKmerNode.put(node.firstKmer(), list);
		}
		list.add(node);
	}
	private void removeFromGraph(KmerPathNodeKmerNode node) {
		Collection<KmerPathNodeKmerNode> list = graphByKmerNode.get(node.firstKmer());
		if (list == null) return;
		list.remove(node);
		if (list.size() == 0) {
			graphByKmerNode.remove(node.firstKmer());
		}
	}
	public boolean sanityCheck() {
		graphByKmerNode.long2ObjectEntrySet().stream().flatMap(e -> e.getValue().stream()).forEach(kn -> { 
			assert(kn.node().isValid());
			assert(graphByPosition.contains(kn.node()));
		});
		for (KmerPathNode n : graphByPosition) {
			assert(n.isValid());
			assert(evidenceTracker.matchesExpected(new KmerPathSubnode(n)));
		}
		if (Defaults.SANITY_CHECK_MEMOIZATION && Defaults.SANITY_CHECK_MEMOIZATION_ALL_OPERATIONS) {
			if (bestContigCaller != null) assert(bestContigCaller.sanityCheck());
		}
		return true;
	}
	public boolean sanityCheckDisjointNodeIntervals() {
		Map<Long, List<KmerPathNode>> byKmer = graphByPosition
	            .stream()
	            .collect(Collectors.groupingBy(KmerPathNode::firstKmer));
		for (List<KmerPathNode> list : byKmer.values()) {
			if (list.size() == 1) continue;
			ArrayList<KmerPathNode> al = Lists.newArrayList(list);
			al.sort(KmerNodeUtil.ByFirstStart);
			for (int i = 1; i < al.size(); i++) {
				assert(al.get(i - 1).firstEnd() < al.get(i).firstStart());
			}
		}
		return true;
	}
	public int tracking_activeNodes() {
		return graphByPosition.size();
	}
	public int tracking_maxKmerActiveNodeCount() {
		return graphByKmerNode.values().stream().mapToInt(x -> x.size()).max().orElse(0);
	}
	public long tracking_underlyingConsumed() {
		return consumed;
	}
	public int tracking_inputPosition() {
		return nextPosition();
	}
	public int tracking_firstPosition() {
		if (graphByPosition.size() == 0) return Integer.MAX_VALUE;
		return graphByPosition.first().firstStart();
	}
	public PositionalDeBruijnGraphTracker getExportTracker() {
		return exportTracker;
	}
	public void setExportTracker(PositionalDeBruijnGraphTracker exportTracker) {
		this.exportTracker = exportTracker;
	}
	public AssemblyChunkTelemetry getTelemetry() {
		return telemetry;
	}
	public void setTelemetry(AssemblyChunkTelemetry telemetry) {
		this.telemetry = telemetry;
	}
	public ContigStats tracking_lastContig() {
		return stats;
	}
	public int tracking_contigsCalled() {
		return contigsCalled;
	}
}
