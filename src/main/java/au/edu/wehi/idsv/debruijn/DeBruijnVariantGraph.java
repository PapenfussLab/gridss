package au.edu.wehi.idsv.debruijn;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Log;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import org.apache.commons.lang3.ArrayUtils;

import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.NonReferenceReadPair;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.RemoteEvidence;
import au.edu.wehi.idsv.SAMEvidenceSource;
import au.edu.wehi.idsv.SoftClipEvidence;

import com.google.common.collect.Sets;

public abstract class DeBruijnVariantGraph<T extends DeBruijnNodeBase> extends DeBruijnGraphBase<T> {
	/**
	 * Once a contig has support from this many distinct evidence sources, the exact values
	 * are unlikely to matter too much and with 16m support kmers in DREAM data, large support
	 * sets are quite computationally expensive
	 */
	private static final int SUPPORT_SIZE_TO_START_APPROXIMATION = 1024;
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
		VariantEvidence graphEvidence = VariantEvidence.createRemoteReadEvidence(getK(), pair);
		addEvidenceKmers(graphEvidence);
	}
	public void addEvidence(SoftClipEvidence read) {
		if (read instanceof RemoteEvidence && !processContext.getAssemblyParameters().includeRemoteSoftClips) {
			// ignore remote soft clips
			return;
		}
		VariantEvidence graphEvidence = VariantEvidence.createSoftClipEvidence(getK(), read);
		addEvidenceKmers(graphEvidence);
	}
	public void removeEvidence(NonReferenceReadPair pair) {
		VariantEvidence graphEvidence = VariantEvidence.createRemoteReadEvidence(getK(), pair);
		removeEvidenceKmers(graphEvidence);
	}
	public void removeEvidence(SoftClipEvidence read) {
		VariantEvidence graphEvidence = VariantEvidence.createSoftClipEvidence(getK(), read);
		removeEvidenceKmers(graphEvidence);
	}
	protected void addEvidenceKmers(VariantEvidence evidence) {
		int readKmerOffset = 0;
		SAMRecord record = evidence.getSAMRecord();
		for (ReadKmer readKmer : new ReadKmerIterable(getK(), record.getReadBases(), record.getBaseQualities(), evidence.isReversed(), evidence.isComplemented())) {
			if (!shouldSkipKmer(evidence, readKmerOffset, readKmer)) {
				T node = createNode(evidence, readKmerOffset, readKmer);
				T graphNode = add(readKmer.kmer, node);
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
		SAMRecord record = evidence.getSAMRecord();
		for (ReadKmer readKmer : new ReadKmerIterable(getK(), record.getReadBases(), record.getBaseQualities(), evidence.isReversed(), evidence.isComplemented())) {
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
	@Override
	public byte[] getBaseCalls(List<Long> path) {
		byte[] bases = super.getBaseCalls(path); 
		return bases;
	}
	@Override
	public byte[] getBaseQuals(List<Long> path) {
		byte[] quals = super.getBaseQuals(path); 
		return quals;
	}
//	/**
//	 * Gets the length of the longest soft clip in the given contig
//	 * @param supportingReads
//	 * @return length of longest read
//	 */
//	protected int getMaxSoftClipLength(List<Long> path, int anchorLength) {
//		int readLength = 0;
//		int offset = 0;
//		for (Long kmer : path) {
//			// don't consider SC length of reference kmers - they may be SC support for a different breakend 
//			if (offset >= anchorLength) {
//				for (SAMRecord r : getKmer(kmer).getSupportingReads()) {
//					readLength = Math.max(readLength, direction == BreakendDirection.Forward ? SAMRecordUtil.getEndSoftClipLength(r) : SAMRecordUtil.getStartSoftClipLength(r));
//				}
//			}
//			offset++;
//		}
//		return readLength;
//	}
	protected class ContigAssembly {
		public Set<DirectedEvidence> support;
		public byte[] baseCalls;
		public byte[] baseQuals;
		public int normalBaseCount;
		public int tumourBaseCount;
		public boolean containsAnchoredSupport() {
			for (DirectedEvidence e : support) {
				if (e instanceof SoftClipEvidence) return true;
			}
			return false;
		}
	}
	protected ContigAssembly debruijnContigAssembly(List<Long> path, int referenceKmersAtStartOfPath) {
		//return assembleMaintainable(path, referenceKmersAtStartOfPath);
		return assembleSinglePass(path, referenceKmersAtStartOfPath);
	}
	@SuppressWarnings("unused")
	private ContigAssembly assembleMaintainable(final List<Long> path, final int referenceKmersAtStartOfPath) {
		final List<Long> breakendPath = path.subList(Math.min(referenceKmersAtStartOfPath, path.size()), path.size());
		return new ContigAssembly() {{
			support = getSupportingEvidence(breakendPath);
			baseCalls = getBaseCalls(path);
			baseQuals = getBaseQuals(path);
			normalBaseCount = getEvidenceKmerCount(breakendPath, false);
			tumourBaseCount = getEvidenceKmerCount(breakendPath, true);
			if (breakendPath.size() > 0) {
				// convert kmer counts into base counts
				// this approximation assumes that read do not deviate then rejoin the assembly path and no misassemblies
				normalBaseCount += (getK() - 1) * (getSupportCount(support, false) - getSupportCount(getKmer(path.get(0)).getSupportingEvidenceList(), false));
				tumourBaseCount += (getK() - 1) * (getSupportCount(support, true) - getSupportCount(getKmer(path.get(0)).getSupportingEvidenceList(), true));
			}
		}};
	}
	private int getSupportCount(Iterable<DirectedEvidence> support, boolean countTumour) {
		int count = 0;
		for (DirectedEvidence e : support) {
			if (((SAMEvidenceSource)e.getEvidenceSource()).isTumour() == countTumour) {
				count++;
			}
		}
		return count;
	}
	private ContigAssembly assembleSinglePass(List<Long> path, int referenceKmersAtStartOfPath, int referenceKmersAtEndOfPath) {
		ContigAssembly ca = new ContigAssembly();
		ca.baseCalls = getBaseCalls(path);
		ca.support = Sets.newHashSet();
		ca.tumourBaseCount = 0;
		ca.normalBaseCount = 0;
		List<Integer> qual = new ArrayList<Integer>(path.size());
		int offset = 0;
		int nAnchorSupportCount = 0;
		int tAnchorSupportCount = 0;
		for (Long node : path) {
			T kmer = getKmer(node);
			List<DirectedEvidence> list = kmer.getSupportingEvidenceList();
			qual.add(kmer.getWeight() - list.size());
			if (offset >= referenceKmersAtStartOfPath) {
				for (DirectedEvidence e : list) {
					SAMEvidenceSource source = (SAMEvidenceSource)e.getEvidenceSource();
					boolean evidenceIsTumour = source != null && source.isTumour();
					if (evidenceIsTumour) {
						ca.tumourBaseCount++;
					} else {
						ca.normalBaseCount++;
					}
				}
				if (ca.support.size() < SUPPORT_SIZE_TO_START_APPROXIMATION) {
					ca.support.addAll(list);
					if (ca.support.size() > SUPPORT_SIZE_HARD_LIMIT) {
						log.warn(String.format("Hit support size hard limit of %d - no longer processing additional support", SUPPORT_SIZE_HARD_LIMIT));
					}
				} else if (ca.support.size() < SUPPORT_SIZE_HARD_LIMIT) {
					// TODO: try a Bloom filter on the support set instead
					// approximate our level of support
					if (offset % (getK() / 2) == 0) {
						// reduce our update frequency to twice per kmer bases
						// this should still get most reads as soft clips
						// are added as the first evidence and RPs should
						// map fully
						ca.support.addAll(list);
						if (ca.support.size() > SUPPORT_SIZE_HARD_LIMIT) {
							log.warn(String.format("Hit support size hard limit of %d - no longer processing additional support", SUPPORT_SIZE_HARD_LIMIT));
						}
					}
				} else {
					// already reached hard limit
				}
				if (offset == referenceKmersAtStartOfPath) {
					nAnchorSupportCount = ca.normalBaseCount;
					tAnchorSupportCount = ca.tumourBaseCount;
				}
			}
			offset++;
		}
		// pad out qualities to match the path length
		for (int i = 0; i < getK() - 1; i++) qual.add(qual.get(qual.size() - 1));
		ca.baseQuals = rescaleBaseQualities(qual);
		// adjust base counts for non-anchored reads
		int nSupport = 0;
		int tSupport = 0;
		for (DirectedEvidence e : ca.support) {
			if (((SAMEvidenceSource)e.getEvidenceSource()).isTumour()) {
				tSupport++;
			} else {
				nSupport++;
			}
		}
		ca.normalBaseCount += (getK() - 1) * (nSupport - nAnchorSupportCount);
		ca.tumourBaseCount += (getK() - 1) * (tSupport - tAnchorSupportCount);
		return ca;
	}
}
