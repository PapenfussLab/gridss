package au.edu.wehi.idsv;

import java.util.Arrays;
import java.util.List;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecord.SAMTagAndValue;
import htsjdk.samtools.util.SequenceUtil;

import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.lang3.tuple.Pair;

import au.edu.wehi.idsv.sam.CigarUtil;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.sam.SamTags;

import com.google.api.client.util.Lists;
import com.google.common.collect.ImmutableList;

public class RealignedRemoteSAMRecordAssemblyEvidence extends RealignedSAMRecordAssemblyEvidence implements RemoteEvidence {
	private final RealignedSAMRecordAssemblyEvidence local;
	public RealignedRemoteSAMRecordAssemblyEvidence(RealignedSAMRecordAssemblyEvidence assembly) {
		this(assembly, remoteRecords(assembly.getEvidenceSource().getContext(), assembly));
	}
	private RealignedRemoteSAMRecordAssemblyEvidence(RealignedSAMRecordAssemblyEvidence assembly, Pair<SAMRecord, SAMRecord> remotes) {
		super(assembly.getEvidenceSource(), remotes.getLeft(), ImmutableList.of(remotes.getRight()));
		this.local = assembly;
	}
	/**
	 * Translate assembly SAMRecord so the realigned record is considered the anchored one
	 * @param assembly
	 * @return
	 */
	private static Pair<SAMRecord, SAMRecord> remoteRecords(ProcessingContext pc, RealignedSAMRecordAssemblyEvidence assembly) {
		SAMRecord realign = SAMRecordUtil.clone(assembly.getRemoteSAMRecord());
		SAMRecord anchor = SAMRecordUtil.clone(assembly.getSAMRecord());
		assert(!realign.getReadUnmappedFlag());
		int realignLength = realign.getReadLength();
		int anchorLength = anchor.getReadLength();
		BreakendDirection direction;
		if ((assembly.getBreakendSummary().direction == BreakendDirection.Forward && realign.getReadNegativeStrandFlag()) || 
				(assembly.getBreakendSummary().direction == BreakendDirection.Backward && !realign.getReadNegativeStrandFlag())) {
			direction = BreakendDirection.Forward;
			realignLength -= SAMRecordUtil.getEndSoftClipLength(realign);
		} else {
			direction = BreakendDirection.Backward;
			realignLength -= SAMRecordUtil.getStartSoftClipLength(realign);
		}
		realign.setAttribute(SamTags.ASSEMBLY_DIRECTION, direction.toChar());
		anchor.setAttribute(SamTags.ASSEMBLY_DIRECTION, null);
		if (realign.getReadNegativeStrandFlag()) {
			// assembly anchors are always +ve strand
			realign.setReadNegativeStrandFlag(!realign.getReadNegativeStrandFlag());
			anchor.setReadNegativeStrandFlag(!anchor.getReadNegativeStrandFlag());
		}
		// swap sequences (AAARR, RR) -> (RRR, RRRAA)
		realign.setReadBases(anchor.getReadBases().clone());
		realign.setBaseQualities(anchor.getBaseQualities().clone());
		if (anchor.getReadNegativeStrandFlag()) {
			// adjust for revcomp 
			SequenceUtil.reverseComplement(realign.getReadBases());
			ArrayUtils.reverse(realign.getBaseQualities());
		}
		// update CIGARS:
		List<CigarElement> cl = Lists.newArrayList(realign.getCigar().getCigarElements());
		cl.add(direction == BreakendDirection.Forward ? cl.size() : 0, new CigarElement(anchorLength - realignLength, CigarOperator.S));
		CigarUtil.clean(cl);
		realign.setCigar(new Cigar(cl));
		if ((direction == BreakendDirection.Forward && !anchor.getReadNegativeStrandFlag()) || 
				direction == BreakendDirection.Backward && anchor.getReadNegativeStrandFlag()) {
			anchor.setCigar(CigarUtil.trimReadBases(anchor.getCigar(), realignLength, 0));
			anchor.setReadBases(Arrays.copyOfRange(anchor.getReadBases(), realignLength, anchor.getReadBases().length));
			anchor.setBaseQualities(Arrays.copyOfRange(anchor.getBaseQualities(), realignLength, anchor.getBaseQualities().length));
		} else {
			anchor.setCigar(CigarUtil.trimReadBases(anchor.getCigar(), 0, realignLength));
			anchor.setReadBases(Arrays.copyOfRange(anchor.getReadBases(), 0, anchor.getReadBases().length - realignLength));
			anchor.setBaseQualities(Arrays.copyOfRange(anchor.getBaseQualities(), 0, anchor.getBaseQualities().length - realignLength));
		}
		
		realign.setFirstOfPairFlag(true);
		anchor.setFirstOfPairFlag(false);
		realign.setSecondOfPairFlag(false);
		anchor.setSecondOfPairFlag(true);
		// Copy assembly annotations across
		for (SAMTagAndValue attribute : anchor.getAttributes()) {
			realign.setAttribute(attribute.tag, attribute.value);
		}
		
		// Sanity checks
		assert(BreakpointFastqEncoding.getEncodedBreakendOffset(anchor.getReadName()) == 0);
		assert(realign.getMappingQuality() >= pc.getConfig().minReadMapq);
		assert(anchor.getMappingQuality() >= pc.getConfig().minReadMapq);
		
		// realign is now the anchor
		return Pair.of(realign, anchor);
	}
	public RealignedSAMRecordAssemblyEvidence asLocal() {
		return local;
	}
	@Override
	public String getEvidenceID() {
		return "R" + super.getEvidenceID();
	}
	@Override
	public String toString() {
		return "R" + super.toString();
	}
}
