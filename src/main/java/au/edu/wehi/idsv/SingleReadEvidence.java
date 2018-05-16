package au.edu.wehi.idsv;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Streams;

import au.edu.wehi.idsv.debruijn.ContigCategorySupportHelper;
import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.sam.SamTags;
import au.edu.wehi.idsv.util.IntervalUtil;
import au.edu.wehi.idsv.util.MessageThrottler;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SequenceUtil;

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
	private final List<Boolean> categorySupport;
	private String evidenceid;
	private boolean unableToCalculateHomology = false;
	private String associatedAssemblyName;
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
			if (!hasForwardSR && SAMRecordUtil.getEndSoftClipLength(record) > 0) {
				list.add(SoftClipEvidence.create(source, BreakendDirection.Forward, record));
			}
			if (!hasBackwardSR && SAMRecordUtil.getStartSoftClipLength(record) > 0) {
				list.add(SoftClipEvidence.create(source, BreakendDirection.Backward, record));
			}
			list.addAll(IndelEvidence.create(source, minIndelSize, record));
		} catch (IllegalArgumentException iae) {
			if (!MessageThrottler.Current.shouldSupress(log, "SingleReadEvidence.createEvidence() failure")) {
				String msg = String.format("createEvidence(): Error processing %s from %s. Ignoring read. "
						+ "This should not happen and may be caused by a internally inconsistent input file (e.g. SA tag does not match read alignment). "
						+ "Please raise an issue at https://github.com/PapenfussLab/gridss/issues "
						+ "and include this stack trace as well as offending SAM record.", record.getReadName(), source == null ? null : source.getFile());
				log.error(iae, msg);
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
		assert(offsetLocalEnd == offsetUnmappedStart || offsetLocalStart == offsetUnmappedEnd);
		assert(offsetUnmappedEnd == offsetRemoteStart || offsetRemoteEnd == offsetUnmappedStart);
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
				assert(offsetLocalEnd - offsetLocalStart == Math.min(2, localInexactMargin));
				if (offsetUnmappedStart >= offsetLocalStart) {
					offsetLocalStart = offsetLocalEnd;
				} else {
					offsetLocalEnd = offsetLocalStart;
				}
			}
			if (remoteInexactMargin > 0) {
				assert(offsetRemoteEnd - offsetRemoteStart == Math.min(2, remoteInexactMargin));
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
		this.categorySupport = getCategorySupportBreakdown(offsetLocalStart, offsetLocalEnd);
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
	@Override
	public boolean isFromMultimappingFragment() {
		return record.getAttribute(SamTags.MULTIMAPPING_FRAGMENT) != null;
	}
	/**
	 * Determines whether this evidence involves the primary read alignment
	 * @return true if the primary read alignment supports this evidence,
	 * false if supported by only supplementary alignments
	 */
	public boolean involvesPrimaryReadAlignment() {
		return !getSAMRecord().getSupplementaryAlignmentFlag();
	}
	/**
	 * Determines which categories 
	 * @param offsetLocalEnd 
	 * @param offsetLocalStart 
	 * @return
	 */
	protected List<Boolean> getCategorySupportBreakdown(int offsetLocalStart, int offsetLocalEnd) {
		if (getEvidenceSource() instanceof AssemblyEvidenceSource) {
			if (!isUnanchored && getSAMRecord().hasAttribute(SamTags.ASSEMBLY_CATEGORY_COVERAGE_CIGAR) && getEvidenceSource().getContext().getConfig().getVariantCalling().requireAssemblyCategorySupport) {
				// Example: alignment is no negative strand
				// 5S 10M 15S
				//   ^ evaluting this breakend
				// position = read base position immediately after breakend
				int nominalPosition = location.direction == BreakendDirection.Backward ? offsetLocalStart : offsetLocalEnd + 1;
				int leftPos = nominalPosition - (location.nominal - location.start);
				int rightPos = nominalPosition + (location.end - location.nominal);
				if (getSAMRecord().getReadNegativeStrandFlag()) {
					// no - 1 since before position -> after position when flipping
					leftPos = getSAMRecord().getReadLength() - leftPos;
					rightPos = getSAMRecord().getReadLength() - rightPos;
				}
				String cigar = getSAMRecord().getStringAttribute(SamTags.ASSEMBLY_CATEGORY_COVERAGE_CIGAR);
				// needs support at both extremes of the homology
				List<Boolean> categories = Streams.zip(
						ContigCategorySupportHelper.supportsBreakendBefore(leftPos, cigar).stream(),
						ContigCategorySupportHelper.supportsBreakendBefore(rightPos, cigar).stream(),
						(a, b) -> (Boolean)(a && b))
					.collect(Collectors.toList());
				// CIGAR only include categories that actually have evidence. Need to pad out to number of categories we actually have 
				while (categories.size() < getEvidenceSource().getContext().getCategoryCount()) {
					categories.add(false);
				}
				return categories;
			} else {
				return Stream.generate(() -> true)
					.limit(getEvidenceSource().getContext().getCategoryCount())
					.collect(Collectors.toList());
			}
		} else {
			if (getEvidenceSource() == null) return null;
			List<Boolean> list = Stream.generate(() -> false)
					.limit(Math.max(getEvidenceSource().getContext().getCategoryCount(), getEvidenceSource().getSourceCategory() + 1))
					.collect(Collectors.toList());
			list.set(getEvidenceSource().getSourceCategory(), true);
			return list;
		}
	}
	public List<Boolean> getCategorySupportBreakdown() {
		return categorySupport;
	}
	@Override
	public List<String> getOriginatingFragmentID(int category) {
		if (AssemblyAttributes.isAssembly(getSAMRecord())) {
			return new AssemblyAttributes(record).getOriginatingFragmentID(category);
		}
		return source.getSourceCategory() == category ? ImmutableList.of(record.getReadName()) : ImmutableList.of();
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
			return new AssemblyAttributes(record).getAssemblyTotalReadSupportCount();
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
}