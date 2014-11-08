package au.edu.wehi.idsv.debruijn;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Log;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.lang3.ArrayUtils;

import au.edu.wehi.idsv.AssemblyBuilder;
import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.NonReferenceReadPair;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.RealignedRemoteSoftClipEvidence;
import au.edu.wehi.idsv.SAMEvidenceSource;
import au.edu.wehi.idsv.SoftClipEvidence;

public abstract class DeBruijnVariantGraph<T extends DeBruijnNodeBase> extends DeBruijnGraphBase<T> {
	/**
	 * Once a contig has support from this many distinct evidence sources, the exact values
	 * are unlikely to matter too much and with 16m support kmers in DREAM data, large support
	 * sets are quite computationally expensive
	 */
	private static final int SUPPORT_SIZE_TO_START_APPROXIMATION = 1024;
	private static final int SUPPORT_SIZE_HARD_LIMIT = 100000;
	private static final Log log = Log.getInstance(DeBruijnVariantGraph.class);
	protected final BreakendDirection direction;
	protected final ProcessingContext processContext;
	protected final AssemblyEvidenceSource source;
	public DeBruijnVariantGraph(ProcessingContext context, AssemblyEvidenceSource source, int k, BreakendDirection direction) {
		super(k);
		this.processContext = context;
		this.source = source;
		this.direction = direction;
	}
	public BreakendDirection getDirection() {
		return direction;
	}
	public void addEvidence(DirectedEvidence evidence) {
		if (evidence instanceof RealignedRemoteSoftClipEvidence) {
			// TODO: should we assemble these? for now, ignore remote soft clips - as they get assembled on the other side
			// TODO: FIXME: add these in
		} else if (evidence instanceof NonReferenceReadPair) {
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
		VariantEvidence graphEvidence = VariantEvidence.createRemoteReadEvidence(direction, getK(), pair);
		addEvidenceKmers(graphEvidence);
	}
	public void addEvidence(SoftClipEvidence read) {
		VariantEvidence graphEvidence = VariantEvidence.createSoftClipEvidence(direction, getK(), read);
		addEvidenceKmers(graphEvidence);
	}
	public void removeEvidence(NonReferenceReadPair pair) {
		VariantEvidence graphEvidence = VariantEvidence.createRemoteReadEvidence(direction, getK(), pair);
		removeEvidenceKmers(graphEvidence);
	}
	public void removeEvidence(SoftClipEvidence read) {
		VariantEvidence graphEvidence = VariantEvidence.createSoftClipEvidence(direction, getK(), read);
		removeEvidenceKmers(graphEvidence);
	}
	protected void addEvidenceKmers(VariantEvidence evidence) {
		int readKmerOffset = 0;
		SAMRecord record = evidence.getSAMRecord();
		for (ReadKmer kmer : new ReadKmerIterable(getK(), record.getReadBases(), record.getBaseQualities(), evidence.isReversed(), evidence.isComplemented())) {
			T node = createNode(evidence, readKmerOffset, kmer);
			if (evidence.isSkippedKmer(readKmerOffset)) {
				// do nothing with skipped kmers
			} else if (kmer.containsAmbiguousBases) {
				// do nothing if the kmer contains an ambiguous base 
			} else {
				T graphNode = add(kmer.kmer, node);
				onEvidenceAdded(graphNode, node, evidence, readKmerOffset, kmer);
			}
			readKmerOffset++;
		}
	}
	protected void removeEvidenceKmers(VariantEvidence evidence) {
		int readKmerOffset = 0;
		SAMRecord record = evidence.getSAMRecord();
		for (ReadKmer kmer : new ReadKmerIterable(getK(), record.getReadBases(), record.getBaseQualities(), evidence.isReversed(), evidence.isComplemented())) {
			T node = createNode(evidence, readKmerOffset, kmer);
			if (evidence.isSkippedKmer(readKmerOffset)) {
				// do nothing with skipped kmers
			} else {
				T graphNode = remove(kmer.kmer, node);
				onEvidenceRemoved(graphNode, node, evidence, readKmerOffset, kmer);
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
		if (direction == BreakendDirection.Backward) {
			ArrayUtils.reverse(bases);
		}
		return bases;
	}
	@Override
	public byte[] getBaseQuals(List<Long> path) {
		byte[] quals = super.getBaseQuals(path); 
		if (direction == BreakendDirection.Backward) {
			ArrayUtils.reverse(quals);
		}
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
	protected AssemblyBuilder debruijnContigAssembly(List<Long> path, int referenceKmersAtStartOfPath) {
		AssemblyBuilder builder = new AssemblyBuilder(processContext, source)
			.direction(direction)
			.assemblyBases(getBaseCalls(path));
		addEvidenceSupportSinglePass(path, referenceKmersAtStartOfPath, builder);
		return builder;
	}
	@SuppressWarnings("unused")
	private AssemblyBuilder addEvidenceSupportMaintainable(List<Long> path, int referenceKmersAtStartOfPath, AssemblyBuilder builder) {
		List<Long> breakendPath = path.subList(referenceKmersAtStartOfPath, path.size());
		Set<DirectedEvidence> support = getSupportingEvidence(breakendPath);
		return builder
				.assemblyBaseQuality(getBaseQuals(path))
				.contributingEvidence(support)
				.assembledBaseCount(getEvidenceBaseCount(breakendPath, false), getEvidenceBaseCount(breakendPath, true));
	}
	private AssemblyBuilder addEvidenceSupportSinglePass(List<Long> path, int referenceKmersAtStartOfPath, AssemblyBuilder builder) {
		Set<DirectedEvidence> support = new HashSet<>();
		int tumourBaseCount = 0;
		int normalBaseCount = 0;
		List<Integer> qual = new ArrayList<Integer>(path.size());
		int offset = 0;
		for (Long node : path) {
			T kmer = getKmer(node);
			List<DirectedEvidence> list = kmer.getSupportingEvidenceList();
			qual.add(kmer.getWeight() - list.size());
			if (offset >= referenceKmersAtStartOfPath) {
				for (DirectedEvidence e : list) {
					SAMEvidenceSource source = (SAMEvidenceSource)e.getEvidenceSource();
					boolean evidenceIsTumour = source != null && source.isTumour();
					if (evidenceIsTumour) {
						tumourBaseCount++;
					} else {
						normalBaseCount++;
					}
				}
				if (support.size() < SUPPORT_SIZE_TO_START_APPROXIMATION) {
					support.addAll(list);
					if (support.size() > SUPPORT_SIZE_HARD_LIMIT) {
						log.warn(String.format("Hit support size hard limit of %d - no longer processing additional support", SUPPORT_SIZE_HARD_LIMIT));
					}
				} else if (support.size() < SUPPORT_SIZE_HARD_LIMIT){
					// approximate our level of support
					if (offset % (getK() / 2) == 0) {
						// reduce our update frequency to twice per kmer bases
						// this should still get most reads as soft clips
						// are added as the first evidence and RPs should
						// map fully
						support.addAll(list);
						if (support.size() > SUPPORT_SIZE_HARD_LIMIT) {
							log.warn(String.format("Hit support size hard limit of %d - no longer processing additional support", SUPPORT_SIZE_HARD_LIMIT));
						}
					}
				} else {
					// reached hard limit - we should log this
				}
			}
			offset++;
		}
		// pad out qualities to match the path length
		for (int i = 0; i < getK() - 1; i++) qual.add(qual.get(qual.size() - 1));
		byte[] quals = rescaleBaseQualities(qual);
		return builder.assemblyBaseQuality(quals)
				.contributingEvidence(support)
				.assembledBaseCount(normalBaseCount, tumourBaseCount);
	}
}
