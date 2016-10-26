package au.edu.wehi.idsv;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import au.edu.wehi.idsv.sam.CigarUtil;

import com.google.common.primitives.Bytes;

public class IndelEvidence implements DirectedBreakpoint {
	private final SAMEvidenceSource source;
	private final SAMRecord record;
	private final BreakpointSummary location;
	private final byte[][] preInsertPostReadBases;
	private final byte[][] preInsertPostQualityScores;
	private final List<CigarElement> indel;
	private IndelEvidence remote;
	private String evidenceid;
	private IndelEvidence(SAMEvidenceSource source, SAMRecord record, List<CigarElement> indel, BreakpointSummary location, byte[][] preInsertPostBases, byte[][] preInsertPostQual) {
		this.source = source;
		this.record = record;
		this.indel = indel;
		this.location = location;
		this.preInsertPostReadBases = preInsertPostBases;
		this.preInsertPostQualityScores = preInsertPostQual;
	}
	public static List<IndelEvidence> create(SAMEvidenceSource source, SAMRecord record) {
		if (record.getReadUnmappedFlag() || record.getCigar() == null) return Collections.emptyList();
		List<IndelEvidence> list = new ArrayList<IndelEvidence>(4);
		List<CigarElement> cl = CigarUtil.decodeNegativeDeletion(record.getCigar().getCigarElements());
		for (int indelStartOffset = 0; indelStartOffset < cl.size(); indelStartOffset++) { // ignore indels at start/end of read
			int indelEndOffset = indelStartOffset; // exclusive end offset
			while (cl.get(indelEndOffset).getOperator().isIndelOrSkippedRegion() && indelEndOffset < cl.size()) {
				indelEndOffset++;
			}
			if (indelStartOffset != indelEndOffset &&
					// Exclude indels at start/end of the read
					indelStartOffset > 0 &&
					indelEndOffset < cl.size()) {
				List<CigarElement> pre = cl.subList(0, indelStartOffset);
				List<CigarElement> indel = cl.subList(indelStartOffset, indelEndOffset);
				List<CigarElement> post = cl.subList(indelEndOffset, cl.size());
				int preReadLength = CigarUtil.readLength(pre);
				int indelReadLength = CigarUtil.readLength(indel);
				int postReadLength = CigarUtil.readLength(post);
				int preRefLength = CigarUtil.referenceLength(pre);
				//int indelRefLength = CigarUtil.referenceLength(indel);
				int postRefLength = CigarUtil.referenceLength(post);
				// TODO: adjust bounds based on sequence homology
				BreakpointSummary location = new BreakpointSummary(
						record.getReferenceIndex(), BreakendDirection.Forward, record.getAlignmentStart() + preRefLength - 1, record.getAlignmentStart() + preRefLength - 1,
						record.getReferenceIndex(), BreakendDirection.Backward, record.getAlignmentEnd() - postRefLength + 1, record.getAlignmentEnd() - postRefLength + 1);
				byte[][] preInsertPostBases = new byte[][] {
						Arrays.copyOf(record.getReadBases(), preReadLength),
						Arrays.copyOfRange(record.getReadBases(), preReadLength, preReadLength + indelReadLength),
						Arrays.copyOfRange(record.getReadBases(), preReadLength + indelReadLength, preReadLength + indelReadLength + postReadLength)};
				byte[][] preInsertPostQual = new byte[][] {
						Arrays.copyOf(record.getBaseQualities(), preReadLength),
						Arrays.copyOfRange(record.getBaseQualities(), preReadLength, preReadLength + indelReadLength),
						Arrays.copyOfRange(record.getBaseQualities(), preReadLength + indelReadLength, preReadLength + indelReadLength + postReadLength)};
				IndelEvidence left = new IndelEvidence(source, record, indel, location, preInsertPostBases, preInsertPostQual);
				IndelEvidence right = new IndelEvidence(source, record, indel, location.remoteBreakpoint(), preInsertPostBases, preInsertPostQual);
				left.remote = right;
				right.remote = left;
				list.add(left);
				list.add(right);
			}
		}
		return list;
	}
	private byte[] local(byte[][] preInsertPost) {
		return preInsertPost[location.direction == BreakendDirection.Forward ? 0 : 2];
	}
	private byte[] indel(byte[][] preInsertPost) {
		return preInsertPost[1];
	}
	private byte[] remote(byte[][] preInsertPost) {
		return preInsertPost[location.direction == BreakendDirection.Forward ? 2 : 0];
	}
	@Override
	public float getBreakendQual() {
		CigarElement e = indel.get(0);
		for (int i = 1; i < indel.size(); i++) {
			if (e.getLength() < indel.get(i).getLength()) {
				e = indel.get(i);
			}
		}
		return (float)getEvidenceSource().getContext().getConfig().getScoring().getModel().scoreIndel(getEvidenceSource().getMetrics(), e.getOperator(), e.getLength(), getLocalMapq()) / 2;
	}
	@Override
	public float getBreakpointQual() {
		return getBreakendQual();
	}

	@Override
	public BreakpointSummary getBreakendSummary() {
		return location;
	}

	@Override
	public byte[] getBreakendSequence() {
		if (location.direction == BreakendDirection.Forward) {
			return Bytes.concat(indel(preInsertPostReadBases), remote(preInsertPostReadBases));
		} else {
			return Bytes.concat(remote(preInsertPostReadBases), indel(preInsertPostReadBases));
		}
	}

	@Override
	public byte[] getBreakendQuality() {
		if (location.direction == BreakendDirection.Forward) {
			return Bytes.concat(indel(preInsertPostQualityScores), remote(preInsertPostQualityScores));
		} else {
			return Bytes.concat(remote(preInsertPostQualityScores), indel(preInsertPostQualityScores));
		}
	}

	@Override
	public byte[] getAnchorSequence() {
		return local(preInsertPostReadBases);
	}

	@Override
	public byte[] getAnchorQuality() {
		return local(preInsertPostQualityScores);
	}

	@Override
	public String getEvidenceID() {
		if (evidenceid == null) {
			StringBuilder sb = new StringBuilder();
			sb.append(record.getReadName());
			if (record.getReadPairedFlag() && record.getFirstOfPairFlag()) {
				sb.append("/1");
			}
			if (record.getReadPairedFlag() && record.getSecondOfPairFlag()) {
				sb.append("/2");
			}
			sb.append("(");
			sb.append(getAnchorSequence().length);
			sb.append(",");
			sb.append(new Cigar(indel));
			sb.append(")");
			evidenceid = sb.toString();
		}
		return evidenceid;
	}

	@Override
	public SAMEvidenceSource getEvidenceSource() {
		return source;
	}

	@Override
	public int getLocalMapq() {
		return record.getMappingQuality();
	}

	@Override
	public int getLocalBaseLength() {
		return local(preInsertPostReadBases).length;
	}

	@Override
	public boolean isBreakendExact() {
		return true;
	}
	@Override
	public int getRemoteMapq() {
		return record.getMappingQuality();
	}
	@Override
	public String getUntemplatedSequence() {
		return new String(indel(preInsertPostReadBases));
	}
	@Override
	public String getHomologySequence() {
		throw new RuntimeException("NYI");
	}
	@Override
	public int getHomologyAnchoredBaseCount() {
		throw new RuntimeException("NYI");
	}
	@Override
	public DirectedBreakpoint asRemote() {
		assert(remote != null);
		return remote;
	}

}
