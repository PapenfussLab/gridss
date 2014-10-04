package au.edu.wehi.idsv.debruijn;

import htsjdk.samtools.SAMRecord;

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
import au.edu.wehi.idsv.SoftClipEvidence;

public abstract class DeBruijnVariantGraph<T extends DeBruijnNodeBase> extends DeBruijnGraphBase<T> {
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
	protected AssemblyBuilder debruijnContigAssembly(List<Long> path) {
		Set<DirectedEvidence> support = getSupportingEvidence(path);
		AssemblyBuilder builder = new AssemblyBuilder(processContext, source)
			.direction(direction)
			.assemblyBases(getBaseCalls(path))
			.assemblyBaseQuality(getBaseQuals(path))
			.contributingEvidence(support)
			.assembledBaseCount(getEvidenceBaseCount(path, false), getEvidenceBaseCount(path, true));
		return builder;
	}
}
