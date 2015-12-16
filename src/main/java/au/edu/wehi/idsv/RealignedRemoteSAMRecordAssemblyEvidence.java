package au.edu.wehi.idsv;

import java.util.List;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
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
		this(assembly, remoteRecords(assembly));
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
	private static Pair<SAMRecord, SAMRecord> remoteRecords(RealignedSAMRecordAssemblyEvidence assembly) {
		SAMRecord realign = SAMRecordUtil.clone(assembly.getRemoteSAMRecord());
		SAMRecord anchor = SAMRecordUtil.clone(assembly.getSAMRecord());
		assert(!realign.getReadUnmappedFlag());
		int realignLength = realign.getReadLength();
		int anchorLength = anchor.getReadLength();
		BreakendDirection direction;
		if ((assembly.getBreakendSummary().direction == BreakendDirection.Forward && realign.getReadNegativeStrandFlag()) || 
				(assembly.getBreakendSummary().direction == BreakendDirection.Backward && !realign.getReadNegativeStrandFlag())) {
			direction = BreakendDirection.Forward;
		} else {
			direction = BreakendDirection.Backward;
		}
		realign.setAttribute(SamTags.ASSEMBLY_DIRECTION, direction.toChar());
		anchor.setAttribute(SamTags.ASSEMBLY_DIRECTION, null);
		if (realign.getReadNegativeStrandFlag()) {
			// assembly anchors are always +ve strand
			realign.setReadNegativeStrandFlag(!realign.getReadNegativeStrandFlag());
			anchor.setReadNegativeStrandFlag(!anchor.getReadNegativeStrandFlag());
		}
		// swap sequences (AAARRR, RRR) -> (RRR, AAARRR)
		byte[] tmp = realign.getReadBases();
		realign.setReadBases(anchor.getReadBases());
		anchor.setReadBases(tmp);
		tmp = realign.getBaseQualities();
		realign.setBaseQualities(anchor.getBaseQualities());
		anchor.setBaseQualities(tmp);
		if (anchor.getReadNegativeStrandFlag()) {
			// adjust for revcomp 
			SequenceUtil.reverseComplement(realign.getReadBases());
			SequenceUtil.reverseComplement(anchor.getReadBases());
			ArrayUtils.reverse(realign.getBaseQualities());
			ArrayUtils.reverse(anchor.getBaseQualities());
		}
		// update CIGARS:
		List<CigarElement> cl = Lists.newArrayList(realign.getCigar().getCigarElements());
		cl.add(direction == BreakendDirection.Forward ? cl.size() : 0, new CigarElement(anchorLength - realignLength, CigarOperator.S));
		CigarUtil.clean(cl);
		realign.setCigar(new Cigar(cl));
		Cigar anchorCigar = anchor.getCigar();
		if ((direction == BreakendDirection.Forward && !anchor.getReadNegativeStrandFlag()) || 
			direction == BreakendDirection.Backward && anchor.getReadNegativeStrandFlag()) {
			anchorCigar = CigarUtil.trimReadBases(anchorCigar, realignLength, 0);
		} else {
			anchorCigar = CigarUtil.trimReadBases(anchorCigar, 0, realignLength);
		}
		anchor.setCigar(anchorCigar);
		
		// TODO: fix read names so realign has an offset (presumably of zero)
		assert(BreakpointFastqEncoding.getEncodedBreakendOffset(anchor.getReadName()) == 0);
		
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
