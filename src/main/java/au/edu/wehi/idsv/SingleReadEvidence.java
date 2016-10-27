package au.edu.wehi.idsv;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.util.IntervalUtil;
import htsjdk.samtools.SAMRecord;

public abstract class SingleReadEvidence implements DirectedEvidence {
	protected static final boolean INCLUDE_CLIPPED_ANCHORING_BASES = false;
	protected final SAMEvidenceSource source;
	private final SAMRecord record;
	private final BreakendSummary location;
	private final int offsetLocalStart;
	private final int offsetLocalEnd;
	private final int offsetUnmappedStart;
	private final int offsetUnmappedEnd;
	private final int offsetRemoteStart;
	private final int offsetRemoteEnd;
	private String evidenceid;
	protected SingleReadEvidence(SAMEvidenceSource source, SAMRecord record, BreakendSummary location,
			int offsetLocalStart, int offsetLocalEnd,
			int offsetUnmappedStart, int offsetUnmappedEnd) {
		this(source, record, location, offsetLocalStart, offsetLocalEnd, offsetUnmappedStart, offsetUnmappedEnd, offsetUnmappedEnd, offsetUnmappedEnd);
	}
	protected SingleReadEvidence(SAMEvidenceSource source, SAMRecord record, BreakendSummary location,
			int offsetLocalStart, int offsetLocalEnd,
			int offsetUnmappedStart, int offsetUnmappedEnd,
			int offsetRemoteStart, int offsetRemoteEnd) {
		// TODO: we could use this to infer the unmapped bounds
		assert(offsetLocalEnd == offsetUnmappedStart || offsetLocalStart == offsetUnmappedEnd);
		assert(offsetUnmappedEnd == offsetUnmappedStart || offsetRemoteEnd == offsetUnmappedStart);
		int remoteBasesToTrim = offsetUnmappedStart - offsetUnmappedEnd;
		if (IntervalUtil.overlapsClosed(offsetLocalStart, offsetLocalEnd - 1, offsetRemoteStart, offsetRemoteEnd - 1)) {
			// Alignments overlap - turns out bwa writes such records
			// To hand these, we'll just reduce the number of bases assigned as remote
			// TODO: update location and consider homologous 
			if (offsetLocalStart < offsetRemoteStart) {
				// local - remote
				offsetUnmappedEnd += remoteBasesToTrim;
				offsetRemoteStart += remoteBasesToTrim;
			} else {
				// remote - local
				offsetRemoteEnd -= remoteBasesToTrim;
				offsetUnmappedStart -= remoteBasesToTrim;
			}
		}
		// validate bounds are valid
		if (offsetLocalStart < 0) throw new IllegalArgumentException();
		if (offsetLocalEnd < 0) throw new IllegalArgumentException();
		if (offsetUnmappedStart < 0) throw new IllegalArgumentException();
		if (offsetLocalEnd < offsetLocalStart) throw new IllegalArgumentException();
		if (offsetUnmappedEnd < offsetUnmappedStart) throw new IllegalArgumentException();
		if (offsetRemoteEnd < offsetRemoteStart) throw new IllegalArgumentException();
		if (offsetLocalEnd > record.getReadLength()) throw new IllegalArgumentException();
		if (offsetUnmappedEnd > record.getReadLength()) throw new IllegalArgumentException();
		if (offsetRemoteEnd > record.getReadLength()) throw new IllegalArgumentException();
		this.source = source;
		this.record = record;
		this.location = location;
		this.offsetLocalStart = offsetLocalStart;
		this.offsetLocalEnd = offsetLocalEnd;
		this.offsetUnmappedStart = offsetUnmappedStart;
		this.offsetUnmappedEnd = offsetUnmappedEnd;
		this.offsetRemoteStart = offsetRemoteStart;
		this.offsetRemoteEnd = offsetRemoteEnd;
		
	}
	
	public SAMRecord getSAMRecord() {
		return record;
	}

	@Override
	public BreakendSummary getBreakendSummary() {
		return location;
	}
	
	private byte[] getBreakend(byte[] b) {
		assert(offsetUnmappedEnd == offsetRemoteStart || offsetUnmappedStart == offsetRemoteEnd);
		return Arrays.copyOfRange(b, Math.min(offsetRemoteStart, offsetUnmappedStart), Math.max(offsetRemoteEnd, offsetUnmappedEnd));
	}

	@Override
	public byte[] getBreakendSequence() {
		return getBreakend(getSAMRecord().getReadBases());
	}

	@Override
	public byte[] getBreakendQuality() {
		return getBreakend(getSAMRecord().getBaseQualities());
	}
	
	private byte[] getAnchor(byte[] b) {
		return Arrays.copyOfRange(b, offsetLocalStart, offsetLocalEnd);
	}

	@Override
	public byte[] getAnchorSequence() {
		return getAnchor(getSAMRecord().getReadBases());
	}

	@Override
	public byte[] getAnchorQuality() {
		return getAnchor(getSAMRecord().getBaseQualities());
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
	public boolean isBreakendExact() {
		return true;
	}

	public String getUntemplatedSequence() {
		return new String(Arrays.copyOfRange(getSAMRecord().getReadBases(), offsetUnmappedStart, offsetUnmappedEnd), StandardCharsets.US_ASCII);
	}
	
	protected void buildEvidenceID(StringBuilder sb) {
		sb.append(record.getReadName());
		if (record.getReadPairedFlag() && record.getFirstOfPairFlag()) {
			sb.append(SAMRecordUtil.FIRST_OF_PAIR_NAME_SUFFIX);
		}
		if (record.getReadPairedFlag() && record.getSecondOfPairFlag()) {
			sb.append(SAMRecordUtil.SECOND_OF_PAIR_NAME_SUFFIX);
		}
	}
	
	@Override
	public String getEvidenceID() {
		if (evidenceid == null) {
			StringBuilder sb = new StringBuilder();
			buildEvidenceID(sb);
			evidenceid = sb.toString();
		}
		return evidenceid;
	}
	
	public String getHomologySequence() {
		throw new RuntimeException("NYI");
	}

	public int getHomologyAnchoredBaseCount() {
		throw new RuntimeException("NYI");
	}
}