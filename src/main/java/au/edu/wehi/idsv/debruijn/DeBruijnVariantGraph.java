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
	protected abstract T createNode(VariantEvidence evidence, int readKmerOffset, ReadKmer kmer);
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
		int readKmerOffset = 0;
		for (ReadKmer readKmer : evidence.getReadKmers()) {
			if (!shouldSkipKmer(evidence, readKmerOffset, readKmer)) {
				T node = createNode(evidence, readKmerOffset, readKmer);
				T graphNode = add(node);
				onEvidenceAdded(graphNode, node, evidence, readKmerOffset, readKmer);
			}
			readKmerOffset++;
		}
	}
	protected boolean shouldSkipKmer(VariantEvidence evidence, int readKmerOffset, ReadKmer readKmer) {
		return evidence.isSkippedKmer(readKmerOffset) || readKmer.containsAmbiguousBases;
	}
	protected void removeEvidenceKmers(VariantEvidence evidence) {
		int readKmerOffset = 0;
		for (ReadKmer readKmer : evidence.getReadKmers()) {
			if (!shouldSkipKmer(evidence, readKmerOffset, readKmer)) {
				T node = createNode(evidence, readKmerOffset, readKmer);
				T graphNode = remove(readKmer.kmer, node);
				onEvidenceRemoved(graphNode, node, evidence, readKmerOffset, readKmer);
			}
			readKmerOffset++;
		}
	}
	/**
	 * Called whenever evidence has been added
	 * @param graphNode de bruijn graph kmer after evidence has been added
	 * @param evidenceNode evidence node that has been added
	 * @param evidence evidence details
	 * @param readKmerOffset kmer offset of this kmer within evidence 
	 * @param kmer evidence kmer
	 */
	protected void onEvidenceAdded(T graphNode, T evidenceNode, VariantEvidence evidence, int readKmerOffset, ReadKmer kmer) { }
	/**
	 * Called whenever evidence has been removed 
	 * @param graphNode de bruijn graph kmer after evidence has been removed
	 * @param evidenceNode evidence that has been removed
	 * @param evidence evidence details
	 * @param readKmerOffset kmer offset of this kmer within evidence 
	 * @param kmer evidence kmer
	 */
	protected void onEvidenceRemoved(T graphNode, T evidenceNode, VariantEvidence evidence, int readKmerOffset, ReadKmer kmer) { }
	private Set<String> getSupport(List<T> path) {
		Set<String> support = Sets.newHashSet();
		// iterate in strides to improve support composition when limit is hit 
		int stride = getK();
		for (int j = 0; j < stride; j++) {
			for (int i = j; i < path.size(); i += stride) {
				T node = path.get(i);
				List<String> nodeSupport = node.getSupportingEvidenceList();
				support.addAll(nodeSupport);
				if (support.size() > SUPPORT_SIZE_HARD_LIMIT) {
					log.warn(String.format("Hit support size hard limit of %d - no longer processing additional support", SUPPORT_SIZE_HARD_LIMIT));
					return support;
				}
			}
		}
		return support;
	}
	private int[] getBaseCountsByCategory(List<T> path, boolean startAnchored, boolean endAnchored) {
		int[] baseCounts = getEvidenceKmerCount(path);
		if (!startAnchored && !endAnchored) {
			// add all the starting bases to the calculation
			int[] firstKmerCounts = path.get(0).getCountByCategory();
			for (int i = 0; i < firstKmerCounts.length; i++) {
				baseCounts[i] += (getK() - 1) * firstKmerCounts[i];
			}
		}
		// TODO: what to do when both anchored? If path.size() < getK() - 1 then expected breakend size is negative;
		// counting kmers instead of bases in this case for now
		return baseCounts;
	}
	protected SAMRecordAssemblyEvidence createAssembly(List<Long> kmers) {
		List<T> path = new ArrayList<T>(kmers.size());
		List<T> breakendPath = new ArrayList<T>(path.size());
		int breakendStartOffset = -1;
		int breakendEndOffset = kmers.size();
		for (int i = 0; i < kmers.size(); i++) {
			T node = getKmer(kmers.get(i));
			path.add(node);
			if (!node.isReference()) {
				breakendPath.add(node);
				if (breakendStartOffset == -1) {
					breakendStartOffset = i;
				}
				breakendEndOffset = i;
			}
		}
		T beforeBreakend = breakendStartOffset > 0 ? path.get(breakendStartOffset -1) : null;
		T afterBreakend = breakendEndOffset + 1 < path.size() ? path.get(breakendEndOffset + 1) : null;
		
		Set<String> breakendSupport = getSupport(breakendPath);
		int[] breakendBaseCounts = getBaseCountsByCategory(breakendPath, beforeBreakend != null, afterBreakend != null);
		// TODO: improve base calling by considering the weight of all kmers contributing to that base position
		byte[] bases = KmerEncodingHelper.baseCalls(kmers, getK());
		byte[] quals = getBaseQuals(kmers);
		SAMRecordAssemblyEvidence ae = null;
		if (beforeBreakend == null && afterBreakend == null) {
			// unanchored
			BreakendSummary breakend = DeBruijnNodeBase.getExpectedBreakend(processContext.getLinear(), breakendPath);
			return AssemblyFactory.createUnanchoredBreakend(processContext, source, breakend, breakendSupport, bases, quals, breakendBaseCounts);
		} else {
			LinearGenomicCoordinate lgc = processContext.getLinear();
			double startBreakendAnchorPosition = DeBruijnNodeBase.getExpectedPositionForDirectAnchor(BreakendDirection.Forward, breakendPath);
			double endBreakendAnchorPosition = DeBruijnNodeBase.getExpectedPositionForDirectAnchor(BreakendDirection.Backward, breakendPath);
			// fall back to expected position of anchor kmer
			// then fall back to expected kmer position - this will give the wrong answer for breakpoint assemblies as it will average the expected location of each breakpoint
			if (Double.isNaN(startBreakendAnchorPosition) && beforeBreakend != null) {
				startBreakendAnchorPosition = beforeBreakend.getExpectedAnchorPosition();
			}
			if (Double.isNaN(endBreakendAnchorPosition) && afterBreakend != null) {
				endBreakendAnchorPosition = afterBreakend.getExpectedAnchorPosition();
			}
			assert(!Double.isNaN(endBreakendAnchorPosition) || !Double.isNaN(startBreakendAnchorPosition));
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
			startBreakendAnchorPosition += getK() - 1;
			
			int startAnchorReferenceIndex = -1;
			int startAnchorPosition =  0;
			int endAnchorReferenceIndex = -1;
			int endAnchorPosition =  0;
			if (!Double.isNaN(startBreakendAnchorPosition)) {
				long pos = (long)Math.round(startBreakendAnchorPosition);
				startAnchorReferenceIndex = lgc.getReferenceIndex(pos); 
				startAnchorPosition = lgc.getReferencePosition(pos);
			}
			if (!Double.isNaN(endBreakendAnchorPosition)) {
				long pos = (long)Math.round(endBreakendAnchorPosition);
				endAnchorReferenceIndex = lgc.getReferenceIndex(pos); 
				endAnchorPosition = lgc.getReferencePosition(pos);
			}
			int startAnchorLength = breakendStartOffset + getK() - 1;
			int endAnchorLength = path.size() + getK() - breakendEndOffset - 2;
			if (endAnchorReferenceIndex == -1) {
				ae = AssemblyFactory.createAnchoredBreakend(processContext, source, BreakendDirection.Forward, breakendSupport,
						startAnchorReferenceIndex, startAnchorPosition, startAnchorLength,
						bases, quals, breakendBaseCounts);
			} else if (startAnchorReferenceIndex == -1) {
				ae = AssemblyFactory.createAnchoredBreakend(processContext, source, BreakendDirection.Backward, breakendSupport,
						endAnchorReferenceIndex, endAnchorPosition, endAnchorLength,
						bases, quals, breakendBaseCounts);
			} else {
				if (startAnchorLength + endAnchorLength > bases.length) {
					// All bases support reference allele - this is likely
					// a) an assembly artifact
					// eg: creates a bubble of nonreference kmers all of which have at least
					// one base supporting the reference
					// MMMMSSS
					// SSSMMMM
					// b) microhomology at breakpoint
					ae = null;
				}
				ae = AssemblyFactory.createAnchoredBreakpoint(processContext, source, breakendSupport,
						startAnchorReferenceIndex, startAnchorPosition, startAnchorLength,
						endAnchorReferenceIndex, endAnchorPosition, endAnchorLength,
						bases, quals, breakendBaseCounts);
			}
		}
		return ae;
	}
}
