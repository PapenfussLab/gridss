package au.edu.wehi.idsv.debruijn;

import htsjdk.samtools.util.Log;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.AssemblyFactory;
import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.BreakendSummary;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.LinearGenomicCoordinate;
import au.edu.wehi.idsv.NonReferenceReadPair;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.RemoteEvidence;
import au.edu.wehi.idsv.SAMRecordAssemblyEvidence;
import au.edu.wehi.idsv.SoftClipEvidence;
import au.edu.wehi.idsv.util.ArrayHelper;

import com.google.common.base.Function;
import com.google.common.collect.Iterables;
import com.google.common.collect.Sets;

public abstract class DeBruijnVariantGraph<T extends DeBruijnNodeBase> extends DeBruijnGraphBase<T> {
	/**
	 * Extremely large support sets are computationally expensive
	 */
	private static final int SUPPORT_SIZE_HARD_LIMIT = 100000;
	private static final Log log = Log.getInstance(DeBruijnVariantGraph.class);
	protected final ProcessingContext processContext;
	protected final AssemblyEvidenceSource source;
	public DeBruijnVariantGraph(ProcessingContext context, AssemblyEvidenceSource source, int k) {
		super(k);
		this.processContext = context;
		this.source = source;
	}
	public void addEvidence(DirectedEvidence evidence) {
		if (evidence instanceof NonReferenceReadPair) {
			addEvidence((NonReferenceReadPair)evidence);
		} else if (evidence instanceof SoftClipEvidence) {
			addEvidence((SoftClipEvidence)evidence);
		} else {
			throw new RuntimeException(String.format("NYI: Unable to add %s evidence to de bruijn graph", evidence));
		}
	}
	protected abstract T createNode(VariantEvidence evidence, int readKmerOffset);
	@Override
	protected T merge(T node, T toAdd) {
		node.add(toAdd);
		return node;
	}
	@Override
	protected T remove(T node, T toRemove) {
		node.remove(toRemove);
		return node;
	}
	public void addEvidence(NonReferenceReadPair pair) {
		VariantEvidence graphEvidence = new VariantEvidence(getK(), pair, processContext.getLinear());
		addEvidenceKmers(graphEvidence);
	}
	public void addEvidence(SoftClipEvidence read) {
		if (read instanceof RemoteEvidence && !processContext.getAssemblyParameters().includeRemoteSoftClips) {
			// ignore remote soft clips
			return;
		}
		VariantEvidence graphEvidence = new VariantEvidence(getK(), read, processContext.getLinear());
		addEvidenceKmers(graphEvidence);
	}
	public void removeEvidence(NonReferenceReadPair pair) {
		VariantEvidence graphEvidence = new VariantEvidence(getK(), pair, processContext.getLinear());
		removeEvidenceKmers(graphEvidence);
	}
	public void removeEvidence(SoftClipEvidence read) {
		VariantEvidence graphEvidence = new VariantEvidence(getK(), read, processContext.getLinear());
		removeEvidenceKmers(graphEvidence);
	}
	protected void addEvidenceKmers(VariantEvidence evidence) {
		PackedReadKmerList kmers = evidence.getKmers();
		for (int i = 0; i < kmers.length(); i++) {
			if (!shouldSkipKmer(evidence, i)) {
				T node = createNode(evidence, i);
				T graphNode = add(node);
				onEvidenceAdded(graphNode, node, evidence, i);
			}
		}
	}
	protected boolean shouldSkipKmer(VariantEvidence evidence, int readKmerOffset) {
		return evidence.isSkippedKmer(readKmerOffset) || evidence.containsAmbiguousBases(readKmerOffset);
	}
	protected void removeEvidenceKmers(VariantEvidence evidence) {
		PackedReadKmerList kmers = evidence.getKmers();
		for (int i = 0; i < kmers.length(); i++) {
			if (!shouldSkipKmer(evidence, i)) {
				T node = createNode(evidence, i);
				T graphNode = remove(kmers.kmer(i), node);
				onEvidenceRemoved(graphNode, node, evidence, i);
			}
		}
	}
	/**
	 * Called whenever evidence has been added
	 * @param graphNode de bruijn graph kmer after evidence has been added
	 * @param evidenceNode evidence node that has been added
	 * @param evidence evidence details
	 * @param readKmerOffset kmer offset of this kmer within evidence 
	 */
	protected void onEvidenceAdded(T graphNode, T evidenceNode, VariantEvidence evidence, int readKmerOffset) { }
	/**
	 * Called whenever evidence has been removed 
	 * @param graphNode de bruijn graph kmer after evidence has been removed
	 * @param evidenceNode evidence that has been removed
	 * @param evidence evidence details
	 * @param readKmerOffset kmer offset of this kmer within evidence 
	 * @param kmer evidence kmer
	 */
	protected void onEvidenceRemoved(T graphNode, T evidenceNode, VariantEvidence evidence, int readKmerOffset) { }
	private Set<String> getSupport(List<List<T>> breakendPathAllKmers) {
		Set<String> support = Sets.newHashSet();
		// iterate in strides to improve support composition when limit is hit 
		int stride = getK();
		for (int j = 0; j < stride; j++) {
			for (int i = j; i < breakendPathAllKmers.size(); i += stride) {
				for (T node : breakendPathAllKmers.get(i)) {
					List<String> nodeSupport = node.getSupportingEvidenceList();
					support.addAll(nodeSupport);
					if (support.size() > SUPPORT_SIZE_HARD_LIMIT) {
						log.warn(String.format("Hit support size hard limit of %d - no longer processing additional support", SUPPORT_SIZE_HARD_LIMIT));
						return support;
					}
				}
			}
		}
		return support;
	}
	private int[] getBaseCountsByCategory(List<List<T>> pathAllKmers, boolean startAnchored, boolean endAnchored) {
		int[] baseCounts = getEvidenceKmerCount(pathAllKmers);
		if (!startAnchored && !endAnchored) {
			// add all the starting bases to the calculation
			int[] firstKmerCounts = new int[0];
			for (T firstKmer : pathAllKmers.get(0)) {
				firstKmerCounts = ArrayHelper.add(firstKmerCounts, firstKmer.getCountByCategory());
			}
			for (int i = 0; i < firstKmerCounts.length; i++) {
				baseCounts[i] += (getK() - 1) * firstKmerCounts[i];
			}
		}
		// TODO: what to do when both anchored? If path.size() < getK() - 1 then expected breakend size is negative;
		// counting kmers instead of bases in this case for now
		return baseCounts;
	}
	protected SAMRecordAssemblyEvidence createAssembly(DeBruijnPathGraph<T, DeBruijnPathNode<T>> pga, List<DeBruijnPathNode<T>> contig) {
		//List<T> path = new ArrayList<T>(contig.size());
		List<List<T>> pathAllKmers = new ArrayList<List<T>>(contig.size());
		//List<T> breakendPath = new ArrayList<T>(contig.size());
		List<List<T>> breakendPathAllKmers = new ArrayList<List<T>>(contig.size());
		for (DeBruijnPathNode<T> pn : contig) {
			boolean isRef = pga.isReference(pn); 
			for (List<T> kmers : pn.getPathAllNodes()) {
				pathAllKmers.add(kmers);
				if (!isRef) {
					breakendPathAllKmers.add(kmers);
				}
			}
		}
		assert(breakendPathAllKmers.size() > 0);
		int breakendStartOffset = pathAllKmers.indexOf(breakendPathAllKmers.get(0));
		int breakendEndOffset = pathAllKmers.indexOf(breakendPathAllKmers.get(breakendPathAllKmers.size() - 1));
		assert(breakendEndOffset - breakendStartOffset == breakendPathAllKmers.size() - 1); // should be a single breakend path
		List<T> beforeBreakend = breakendStartOffset > 0 ? pathAllKmers.get(breakendStartOffset - 1) : null;
		List<T> afterBreakend = breakendEndOffset + 1 < pathAllKmers.size() ? pathAllKmers.get(breakendEndOffset + 1) : null;
		
		if (breakendPathAllKmers.size() < getK() - 1 && beforeBreakend != null && afterBreakend != null) {
			// our bubble is smaller than expected - this is likely due to either:
			// a) assembly artifact: bubble of non-reference kmers all of which have at least one base supporting the reference
			//     MMMMSSS
			//     SSSMMMM			
			// b) microhomology at breakpoint: adjust bases called so we call the middle of the microhomology
			//     making sure we don't overrun either anchor and turn the breakpoint into a breakend
			int basesToAdjust = getK() - 1 - breakendPathAllKmers.size();
			while (basesToAdjust > 0 && (breakendStartOffset > 1 || breakendEndOffset < pathAllKmers.size() - 2)) {
				if (basesToAdjust > 0 && breakendStartOffset > 1) {
					breakendPathAllKmers.add(0, beforeBreakend);
					breakendStartOffset--;
					basesToAdjust--;
					beforeBreakend = pathAllKmers.get(breakendStartOffset - 1);
				}
				if (basesToAdjust > 0 && breakendEndOffset < pathAllKmers.size() - 2) {
					breakendPathAllKmers.add(afterBreakend);
					breakendEndOffset++;
					basesToAdjust--;
					afterBreakend = pathAllKmers.get(breakendEndOffset + 1);
				}
			}
			if (basesToAdjust > 0) {
				// not enough anchored bases either side.
				// TODO: chew into the anchoring kmers and calling with less than k anchor bases either side instead of outright rejection
				log.debug(String.format("Rejecting undersized breakpoint assembly near (%d).", beforeBreakend.get(0).getExpectedPosition()));
				return null;
			}
		}
		Set<String> breakendSupport = getSupport(breakendPathAllKmers);
		int[] breakendBaseCounts = getBaseCountsByCategory(breakendPathAllKmers, beforeBreakend != null, afterBreakend != null);
		byte[] bases = KmerEncodingHelper.baseCalls(KmerEncodingHelper.asKmers(this, Iterables.transform(pathAllKmers, new Function<List<T>, T>() {
			// TODO: improve base calling by considering all kmers contributing to each base position
			@Override
			public T apply(List<T> input) {
				return input.get(0);
			}})), getK());
		byte[] quals = getBaseQuals(pathAllKmers);
		SAMRecordAssemblyEvidence ae = null;
		if (beforeBreakend == null && afterBreakend == null) {
			// unanchored
			BreakendSummary breakend = DeBruijnNodeBase.getExpectedBreakend(processContext.getLinear(), Iterables.concat(breakendPathAllKmers));
			return AssemblyFactory.createUnanchoredBreakend(processContext, source, breakend, breakendSupport, bases, quals, breakendBaseCounts);
		} else {
			LinearGenomicCoordinate lgc = processContext.getLinear();
			Long startBreakendAnchorPosition = null;
			Long endBreakendAnchorPosition = null;
			if (beforeBreakend != null) {
				startBreakendAnchorPosition = DeBruijnNodeBase.getExpectedReferencePosition(beforeBreakend);
			}
			if (afterBreakend != null) {
				endBreakendAnchorPosition = DeBruijnNodeBase.getExpectedReferencePosition(afterBreakend);
			}
			assert(endBreakendAnchorPosition != null || startBreakendAnchorPosition != null);
			// k=3
			// 1234567890
			// MMMSSSMMM
			// 000   | |
			// |111  | |
			// | 222 | |
			// |  333| |
			// |   444 |
			// |    555|
			// |     666
			// ^ start kmer pos
			//       ^ end kmer pos
			// adjust position from start of kmer over to closest reference anchor base position
			if (startBreakendAnchorPosition != null) {
				startBreakendAnchorPosition += getK() - 1;
			}
			
			int startAnchorReferenceIndex = -1;
			int startAnchorPosition =  0;
			int endAnchorReferenceIndex = -1;
			int endAnchorPosition =  0;
			if (startBreakendAnchorPosition != null) {
				startAnchorReferenceIndex = lgc.getReferenceIndex(startBreakendAnchorPosition); 
				startAnchorPosition = lgc.getReferencePosition(startBreakendAnchorPosition);
			}
			if (endBreakendAnchorPosition != null) {
				endAnchorReferenceIndex = lgc.getReferenceIndex(endBreakendAnchorPosition); 
				endAnchorPosition = lgc.getReferencePosition(endBreakendAnchorPosition);
			}
			int startAnchorLength = breakendStartOffset + getK() - 1;
			int endAnchorLength = pathAllKmers.size() + getK() - breakendEndOffset - 2;
			if (endAnchorReferenceIndex == -1) {
				ae = AssemblyFactory.createAnchoredBreakend(processContext, source, BreakendDirection.Forward, breakendSupport,
						startAnchorReferenceIndex, startAnchorPosition, startAnchorLength,
						bases, quals, breakendBaseCounts);
			} else if (startAnchorReferenceIndex == -1) {
				ae = AssemblyFactory.createAnchoredBreakend(processContext, source, BreakendDirection.Backward, breakendSupport,
						endAnchorReferenceIndex, endAnchorPosition, endAnchorLength,
						bases, quals, breakendBaseCounts);
			} else {
				ae = AssemblyFactory.createAnchoredBreakpoint(processContext, source, breakendSupport,
						startAnchorReferenceIndex, startAnchorPosition, startAnchorLength,
						endAnchorReferenceIndex, endAnchorPosition, endAnchorLength,
						bases, quals, breakendBaseCounts);
			}
		}
		return ae;
	}
}
