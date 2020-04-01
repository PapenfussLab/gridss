package au.edu.wehi.idsv;

import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.sam.ChimericAlignment;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.util.IntervalUtil;
import au.edu.wehi.idsv.util.MessageThrottler;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Range;
import htsjdk.samtools.*;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SequenceUtil;

import java.nio.charset.StandardCharsets;
import java.util.*;

public abstract class SingleReadEvidence implements DirectedEvidence {
	private static final Log log = Log.getInstance(SingleReadEvidence.class);
	protected static final boolean INCLUDE_CLIPPED_ANCHORING_BASES = false;
	protected final SAMEvidenceSource source;
	private final SAMRecord record;
	private final BreakendSummary location;
	private final byte[] anchorBases;
	private final byte[] anchorQuals;
	private final String untemplated;
	private final byte[] breakendBases;
	private final byte[] breakendQuals;
	private final boolean isUnanchored;
	/**
	 * Offset in the read alignment of the nominal anchoring base flanking the breakend.
	 * That is, the aligned base immediately next to the clip/indel for which this is evidence for
	 */
	private final int nominalOffset;
	private String evidenceid;
	private boolean unableToCalculateHomology = false;
	private String associatedAssemblyName;
	private int assemblyOffset = Integer.MIN_VALUE;

	public static List<SingleReadEvidence> createEvidence(SAMEvidenceSource source, int minIndelSize, SAMRecord record) {
		if (record.getReadUnmappedFlag()) return Collections.emptyList();
		List<SingleReadEvidence> list = new ArrayList<>(4);
		try {
			List<SplitReadEvidence> srlist = SplitReadEvidence.create(source, record);
			list.addAll(srlist);
			// only add soft clip if there isn't a split read
			boolean hasForwardSR = false;
			boolean hasBackwardSR = false;
			for (SplitReadEvidence sre : srlist) {
				switch (sre.getBreakendSummary().direction) {
					case Forward:
						hasForwardSR = true;
						break;
					case Backward:
						hasBackwardSR = true;
						break;
				}
			}
			// TODO: min clip length
			// not dovetailing
			if (source == null || !SAMRecordUtil.isDovetailing(record, SamPairUtil.PairOrientation.FR, source.getContext().getConfig().dovetailMargin)) {
				if (!hasForwardSR && SAMRecordUtil.getEndSoftClipLength(record) > 0) {

					list.add(SoftClipEvidence.create(source, BreakendDirection.Forward, record));
				}
				if (!hasBackwardSR && SAMRecordUtil.getStartSoftClipLength(record) > 0) {
					list.add(SoftClipEvidence.create(source, BreakendDirection.Backward, record));
				}
			}
			list.addAll(IndelEvidence.create(source, minIndelSize, record));
		} catch (IllegalArgumentException iae) {
			if (!MessageThrottler.Current.shouldSupress(log, "SingleReadEvidence.createEvidence() failure")) {
				boolean sa_error_message_written = false;
				int readLength = record.getReadLength();
				for (ChimericAlignment ca : ChimericAlignment.getChimericAlignments(record)) {
					if (!sa_error_message_written) {
						int len = 0;
						for (CigarElement ce : ca.cigar) {
							if (ce.getOperator().consumesReadBases() || ce.getOperator() == CigarOperator.H) {
								len += ce.getLength();
							}
						}
						if (len != readLength) {
							String msg = String.format("Data sanity check failure: split read alignment of %s from %s have different read lengths. "
									+ " This is typically caused by GATK indel realignment stripping hard clipping from read alignments. "
									+ " Ignoring read", record.getReadName(), source == null ? null : source.getFile());
							log.error(msg);
							sa_error_message_written = true;
						}
					}
				}
				if (!sa_error_message_written) {
					String msg = String.format("createEvidence(): Error processing %s from %s. Ignoring read. "
							+ "This should not happen and may be caused by a internally inconsistent input file (e.g. SA tag does not match read alignment). "
							+ "Please raise an issue at https://github.com/PapenfussLab/gridss/issues "
							+ "and include this stack trace as well as offending SAM record.",
                            record.getReadName(), source == null ? null : source.getFile());
					log.error(iae, msg);
				}
			}
		}
		return list;
	}

	protected SingleReadEvidence(SAMEvidenceSource source, SAMRecord record, BreakendSummary location,
								 int offsetLocalStart, int offsetLocalEnd,
								 int offsetUnmappedStart, int offsetUnmappedEnd,
								 int localInexactMargin) {
		this(source, record, location,
				offsetLocalStart, offsetLocalEnd,
				offsetUnmappedStart, offsetUnmappedEnd,
				offsetLocalStart < offsetUnmappedStart ? offsetUnmappedEnd : offsetUnmappedStart,
				offsetLocalStart < offsetUnmappedStart ? offsetUnmappedEnd : offsetUnmappedStart,
				localInexactMargin, 0);
	}

	protected SingleReadEvidence(SAMEvidenceSource source, SAMRecord record, BreakendSummary location,
								 int offsetLocalStart, int offsetLocalEnd,
								 int offsetUnmappedStart, int offsetUnmappedEnd,
								 int offsetRemoteStart, int offsetRemoteEnd,
								 int localInexactMargin, int remoteInexactMargin) {

		// TODO: we could use this to infer the unmapped bounds
		assert (offsetLocalEnd == offsetUnmappedStart || offsetLocalStart == offsetUnmappedEnd);
		assert (offsetUnmappedEnd == offsetRemoteStart || offsetRemoteEnd == offsetUnmappedStart);
		if (IntervalUtil.overlapsClosed(offsetLocalStart, offsetLocalEnd - 1, offsetRemoteStart, offsetRemoteEnd - 1)) {
			int remoteBasesToTrim = offsetUnmappedStart - offsetUnmappedEnd;
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
		if (localInexactMargin > 0 || remoteInexactMargin > 0) {
			// strip out placeholder anchorings
			// 1X*N1X format has at most 2 anchoring bases
			if (localInexactMargin > 0) {
				assert (offsetLocalEnd - offsetLocalStart == Math.min(2, localInexactMargin));
				if (offsetUnmappedStart >= offsetLocalStart) {
					offsetLocalStart = offsetLocalEnd;
				} else {
					offsetLocalEnd = offsetLocalStart;
				}
			}
			if (remoteInexactMargin > 0) {
				assert (offsetRemoteEnd - offsetRemoteStart == Math.min(2, remoteInexactMargin));
				if (offsetUnmappedEnd >= offsetRemoteEnd) {
					offsetRemoteStart = offsetRemoteEnd;
				} else {
					offsetRemoteEnd = offsetRemoteStart;
				}
			}
			// adjust breakend bounds
			location = location.adjustPosition(Math.max(0, localInexactMargin - 1), Math.max(0, remoteInexactMargin - 1), true);
		} else {
			// TODO calculate microhomology length
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
		if (offsetUnmappedEnd != offsetRemoteStart && offsetUnmappedStart != offsetRemoteEnd) throw new IllegalArgumentException();
		this.source = source;
		this.record = record;
		this.untemplated = new String(Arrays.copyOfRange(record.getReadBases(), offsetUnmappedStart, offsetUnmappedEnd), StandardCharsets.US_ASCII);
		this.anchorBases = Arrays.copyOfRange(record.getReadBases(), offsetLocalStart, offsetLocalEnd);
		this.breakendBases = Arrays.copyOfRange(record.getReadBases(), Math.min(offsetRemoteStart, offsetUnmappedStart), Math.max(offsetRemoteEnd, offsetUnmappedEnd));
		if (record.getBaseQualities() != SAMRecord.NULL_QUALS && record.getBaseQualities() != null) {
			this.anchorQuals = Arrays.copyOfRange(record.getBaseQualities(), offsetLocalStart, offsetLocalEnd);
			this.breakendQuals = Arrays.copyOfRange(record.getBaseQualities(), Math.min(offsetRemoteStart, offsetUnmappedStart), Math.max(offsetRemoteEnd, offsetUnmappedEnd));
		} else {
			this.anchorQuals = null;
			this.breakendQuals = null;
		}
		this.isUnanchored = localInexactMargin > 0 || remoteInexactMargin > 0;
		location = withExactHomology(location);
		if (source != null && source.getContext() != null && source.getContext().getReference() != null && source.getContext().getReference().getSequenceDictionary() != null) {
			location = location.asValidFor(source.getContext().getReference().getSequenceDictionary());
		}
		this.location = location;
		this.nominalOffset = location.direction == BreakendDirection.Forward ? offsetLocalEnd - 1: offsetLocalStart;
	}
	private BreakendSummary withExactHomology(BreakendSummary location) {
		if (!isBreakendExact()) return location;
		// if there is an overlapping split read alignment (such as produced by bwa)
		// then we'll trust the overlap given by the aligner as our homology size
		// even though this could be an inexact homolog
		if (location.end - location.start != 0) return location;
		// If there's inserted sequence that hasn't been aligned to either side then we don't have a homology.
		// Edge case: technically this isn't correct. Sequences such as
		// tandem repeats can have both inserted sequence and sequence homology.
		if (untemplated.length() > 0) return location;
		if (location instanceof BreakpointSummary) {
			if (source != null && source.getContext() != null && source.getContext().getReference() !=  null) {
				ReferenceLookup lookup = source.getContext().getReference();
				BreakpointSummary bp = (BreakpointSummary) location;
				int localBasesMatchingRemoteReference;
				int remoteBasesMatchingLocalReference;
				if (bp.direction == BreakendDirection.Forward) {
					// anchor -> breakend
					remoteBasesMatchingLocalReference = homologyLength(lookup, bp.referenceIndex, bp.nominal + 1, 1, breakendBases, 0, 1);
					if (bp.direction2  == BreakendDirection.Backward) {
						localBasesMatchingRemoteReference = homologyLength(lookup, bp.referenceIndex2, bp.nominal2 - 1, -1, anchorBases, anchorBases.length - 1, -1);
					} else {
						localBasesMatchingRemoteReference = homologyLength(lookup, bp.referenceIndex2, bp.nominal2 + 1,  1, anchorBases, anchorBases.length - 1, -1);
					}
				} else {
					remoteBasesMatchingLocalReference = homologyLength(lookup, bp.referenceIndex, bp.nominal - 1, -1, breakendBases, breakendBases.length - 1, -1);
					if (bp.direction2  == BreakendDirection.Forward) {
						localBasesMatchingRemoteReference = homologyLength(lookup, bp.referenceIndex2, bp.nominal2 + 1,  1, anchorBases, 0, 1);
					} else {
						localBasesMatchingRemoteReference = homologyLength(lookup, bp.referenceIndex2, bp.nominal2 - 1, -1, anchorBases, 0, 1);
					}
				}
				BreakpointSummary adjusted = bp.adjustPosition(localBasesMatchingRemoteReference, remoteBasesMatchingLocalReference, false);
				// push the bounds back in if we overrun our contigs bounds to the homology 
				if (adjusted.start <= 0 && localBasesMatchingRemoteReference > 0) localBasesMatchingRemoteReference--;
				if (adjusted.start2 <= 0 && remoteBasesMatchingLocalReference > 0) remoteBasesMatchingLocalReference--;
				if (adjusted.end > lookup.getSequenceDictionary().getSequence(adjusted.referenceIndex).getSequenceLength() && localBasesMatchingRemoteReference > 0) localBasesMatchingRemoteReference--;
				if (adjusted.end2 > lookup.getSequenceDictionary().getSequence(adjusted.referenceIndex2).getSequenceLength() && remoteBasesMatchingLocalReference > 0) remoteBasesMatchingLocalReference--;
				BreakpointSummary readjusted = bp.adjustPosition(localBasesMatchingRemoteReference, remoteBasesMatchingLocalReference, false);
				return readjusted;
			} else {
				unableToCalculateHomology = true;
				return location;
			}
		} else {
			// not considering unclipping bases since that only affects assembly
			// and is already covered by SAMRecordUtil.unclipExactReferenceMatches()
			return location;
		}
	}
	private static int homologyLength(ReferenceLookup lookup, int referenceIndex, int referencePosition, int referenceStep, byte[] seq, int seqPosition, int seqStep) {
		SAMSequenceRecord refSeq = lookup.getSequenceDictionary().getSequence(referenceIndex);
		int homlen = 0;
		boolean complement = referenceStep != seqStep;
		while (seqPosition >= 0 && seqPosition < seq.length &&
				// next step must still be on the contig
				referencePosition >= 1 && referencePosition <= refSeq.getSequenceLength()) {
			byte base = seq[seqPosition];
			if (complement) {
				base = SequenceUtil.complement(base);
			}
			if (SequenceUtil.basesEqual(base, lookup.getBase(referenceIndex, referencePosition)) &&
					!SequenceUtil.basesEqual(base, SequenceUtil.N)) {
				referencePosition += referenceStep;
				seqPosition += seqStep;
				homlen++;
			} else {
				break;
			}
		}
		return homlen;
	}
	
	public SAMRecord getSAMRecord() {
		return record;
	}

	@Override
	public BreakendSummary getBreakendSummary() {
		return location;
	}
	
	@Override
	public byte[] getBreakendSequence() {
		return breakendBases;
	}

	@Override
	public byte[] getBreakendQuality() {
		return breakendQuals;
	}
	

	@Override
	public byte[] getAnchorSequence() {
		return anchorBases;
	}

	@Override
	public byte[] getAnchorQuality() {
		return anchorQuals;
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
		return !isUnanchored;
	}

	public String getUntemplatedSequence() {
		return untemplated;
	}
	
	protected abstract String getUncachedEvidenceID();
	
	@Override
	public String getEvidenceID() {
		if (evidenceid == null) {
			evidenceid = getUncachedEvidenceID();
		}
		return evidenceid;
	}
	
	public String getHomologySequence() {
		if (unableToCalculateHomology) throw new IllegalStateException("Unable to calculate homology as reference genome has not been supplied");
		if (!isBreakendExact()) return "";
		int homlen = location.end - location.start;
		int locallen = getHomologyAnchoredBaseCount();
		int remotelen = homlen - locallen;
		String strAnchor = new String(anchorBases, StandardCharsets.US_ASCII);
		String strBreakend = new String(breakendBases, StandardCharsets.US_ASCII);
		try {
			if (location.direction == BreakendDirection.Forward) {
				// end of anchor + start of breakend
				return strAnchor.substring(strAnchor.length() - locallen) + strBreakend.substring(0, remotelen);
			} else {
				// end of breakend + start of anchor
				return strBreakend.substring(strBreakend.length() - remotelen) + strAnchor.substring(0, locallen);
			}
		} catch (StringIndexOutOfBoundsException e) {
			String msg = String.format("Sanity check failure: getHomologySequence() failed for %s at %s:%d (%s). Local/remote homology lengths of %d/%d"
					+ " not compatible anchor and breakend lengths of %d/%d",
					record.getReadName(),
					record.getContig(),
					record.getAlignmentStart(),
					getBreakendSummary().toString(source.getContext()),
					locallen,
					remotelen,
					strAnchor.length(),
					strBreakend.length());
			if (!MessageThrottler.Current.shouldSupress(log, "getHomologySequence() failures")) {
				log.error(msg);
			}
			return "";
		}
	}

	public int getHomologyAnchoredBaseCount() {
		if (unableToCalculateHomology) throw new IllegalStateException("Unable to calculate homology as reference genome has not been supplied");
		if (isBreakendExact()) {
			if (location.direction == BreakendDirection.Forward) { 
				return location.nominal - location.start;
			} else  {
				return location.end - location.nominal;
			}
		}
		return 0;
	}
	
	/**
	 * Evidence provides support for no structural variant call
	 * @return true if this evidence is consistent with the reference allele
	 */
	public abstract boolean isReference();
	
	@Override
	public String toString() {
		return getEvidenceID();
	}
	/**
	 * Determines whether this evidence involves the primary read alignment
	 * @return true if the primary read alignment supports this evidence,
	 * false if supported by only supplementary alignments
	 */
	public boolean involvesPrimaryReadAlignment() {
		return !getSAMRecord().getSupplementaryAlignmentFlag();
	}
	@Override
	public Collection<String> getOriginatingFragmentID(int category) {
		if (AssemblyAttributes.isAssembly(getSAMRecord())) {
			return new AssemblyAttributes(record).getOriginatingFragmentID(null, ImmutableSet.of(category), null);
		}
		return source.getSourceCategory() == category ? ImmutableSet.of(record.getReadName()) : ImmutableSet.of();
	}
	@Override
	public double getStrandBias() {
		double bias = 1;
		if (AssemblyAttributes.isAssembly(getSAMRecord())) {
			bias = new AssemblyAttributes(record).getStrandBias();
		}
		if (record.getReadNegativeStrandFlag()) {
			bias = 1 - bias;
		}
		return bias;
	}
	public int constituentReads() {
		if (AssemblyAttributes.isAssembly(getSAMRecord())) {
			AssemblyAttributes aa = new AssemblyAttributes(record);
			return aa.getSupportingReadCount(getBreakendAssemblyContigOffset(), null, null);
		}
		return 1;
	}
	@Override
	public String getAssociatedAssemblyName() {
		if (AssemblyAttributes.isAssembly(getSAMRecord())) {
			return getSAMRecord().getReadName();
		}
		return associatedAssemblyName;
	}
	public void setAssociatedAssemblyName(String associatedAssemblyName) {
		this.associatedAssemblyName = associatedAssemblyName; 
	}

	/**
	 * Position of the breakend relative to the start of the read sequence.
	 * A range is required to account for breakpoint sequence homology.
	 * @return 0-based position relative to the start of the read.
	 * The breakpoint is considered to be immediately prior (in read-space coordinates) to the given position.
	 */
	public Range<Integer> getBreakendAssemblyContigBreakpointInterval() {
		if (isUnanchored) {
			int anchorBases = AssemblyAttributes.getUnanchoredPlacholderAnchoredBases(record);
			if (new AssemblyAttributes(record).getAssemblyDirection() == BreakendDirection.Forward) {
				return Range.closed(anchorBases, anchorBases);
			} else {
				int rl = SAMRecordUtil.getReadLengthIncludingHardClipping(record);
				return Range.closed(rl - anchorBases, rl - anchorBases);
			}
		}
        Range<Integer> r;
		int startAnchoringOffset = nominalOffset + location.start - location.nominal;
		int endAnchoringOffset = nominalOffset + location.end - location.nominal;
		int insertLength = getUntemplatedSequence().length();
		if (location instanceof BreakpointSummary && insertLength > 0) {
			// Inserted sequence - extend the interval over the inserted sequence
			// so we have symmetry across the inserted sequence
			if (location.direction == BreakendDirection.Forward) {
				endAnchoringOffset += insertLength;
			} else {
				startAnchoringOffset -= insertLength;
			}
		}
		if (location.direction == BreakendDirection.Forward) {
			// translate from anchor position coordinates to immediately-before coordinates
			startAnchoringOffset++;
			endAnchoringOffset++;
		}
		if (record.getReadNegativeStrandFlag()) {
			// swap to read-based coordinate system
			r = Range.closed(
					record.getReadLength() - endAnchoringOffset,
					record.getReadLength() -startAnchoringOffset);
		} else {
			r = Range.closed(startAnchoringOffset, endAnchoringOffset);
		}
		return r;
		/*
		boolean onNegative = record.getReadNegativeStrandFlag();
		int nominalBreakendAfterReadOffset = onNegative ? record.getReadLength() - nominalOffset : nominalOffset;
		if (!(location instanceof BreakpointSummary)) {
			// breakend is just the actual position
			assert(location.start == location.nominal);
			assert(location.end == location.nominal);
			r = Range.closed(nominalBreakendAfterReadOffset, nominalBreakendAfterReadOffset);
			return r;
		}
		int homologyBasesAfterNominal =
		// Variant has either inserted sequence or homology
		int insertLength = getUntemplatedSequence().length();
		if (insertLength > 0) {
			if ((location.direction == BreakendDirection.Forward && !record.getReadNegativeStrandFlag()) |
					(location.direction == BreakendDirection.Backward && record.getReadNegativeStrandFlag())) {
				r = Range.closed(nominalBreakendAfterReadOffset, nominalBreakendAfterReadOffset + insertLength);
			} else {
				r = Range.closed(nominalBreakendAfterReadOffset - insertLength + 1, nominalBreakendAfterReadOffset + 1);
			}
		} else {
			if (!record.getReadNegativeStrandFlag()) {
				r = Range.closed(nominalBreakendAfterReadOffset + (location.start - location.nominal), nominalBreakendAfterReadOffset + (location.end - location.nominal));
			} else {
				r = Range.closed(nominalBreakendAfterReadOffset - (location.end - location.nominal), nominalBreakendAfterReadOffset - (location.start - location.nominal));
			}
		}
		return r;
		*/
	}
	public int getBreakendAssemblyContigOffset() {
		if (assemblyOffset == Integer.MIN_VALUE && AssemblyAttributes.isAssembly(record)) {
			AssemblyAttributes aa = new AssemblyAttributes(record);
			assemblyOffset = aa.getMinQualPosition(getBreakendAssemblyContigBreakpointInterval(), null, null);
		}
		return assemblyOffset;
	}
	/**
	 * Offset of the read alignment
	 * @return
	 */
	public int getLocalChimericAlignmentReadOffset() {
		return new ChimericAlignment(record).getFirstAlignedBaseReadOffset();
	}

	@Override
	public SAMRecord getUnderlyingSAMRecord() {
		return getSAMRecord();
	}
}