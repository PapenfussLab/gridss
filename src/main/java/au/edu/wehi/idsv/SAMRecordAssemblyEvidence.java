package au.edu.wehi.idsv;

import gnu.trove.set.hash.TIntHashSet;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.util.Log;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.apache.commons.lang3.StringUtils;

import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.sam.SamTags;
import au.edu.wehi.idsv.vcf.VcfAttributes;
import au.edu.wehi.idsv.vcf.VcfFilter;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.collect.Ordering;
import com.google.common.hash.HashFunction;
import com.google.common.hash.Hashing;
import com.google.common.primitives.Bytes;
import com.google.common.primitives.UnsignedBytes;

public class SAMRecordAssemblyEvidence implements AssemblyEvidence {
	private static final Log log = Log.getInstance(SAMRecordAssemblyEvidence.class);
	private final SAMRecord record;
	private final SAMRecord realignment;
	private final AssemblyEvidenceSource source;
	private final BreakendSummary breakend;
	private final boolean isExact;
	private Collection<DirectedEvidence> evidence = new ArrayList<DirectedEvidence>();
	public SAMRecordAssemblyEvidence(AssemblyEvidenceSource source, SAMRecord assembly, SAMRecord realignment) {
		this.source = source;
		this.record = assembly;
		this.breakend = calculateBreakendFromAlignmentCigar(assembly.getReferenceIndex(), assembly.getAlignmentStart(), assembly.getCigar());
		this.isExact = calculateIsBreakendExactFromCigar(assembly.getCigar());
		this.realignment = realignment == null ? getPlaceholderRealignment() : realignment;
		fixReadPair();
	}
	/**
	 * Lazily calculated evidence IDs of the evidence contributing to this breakend
	 * 
	 * Hashes of the IDs are stored due to the overhead of storing the full strings
	 * and the low likelihood and impact of hash collisions
	 */
	private TIntHashSet evidenceIdhashSet = null;
	private int evidenceSetSize = -1;
	private static HashFunction evidenceIdHash = Hashing.murmur3_32();
	private void ensureEvidenceSet() {
		if (evidenceIdhashSet == null) {
			Object value = record.getAttribute(SamTags.ASSEMLBY_EVIDENCE_HASH);
			if (value instanceof int[]) {
				int[] hashes = (int[])value;
				evidenceIdhashSet = new TIntHashSet(hashes);
				if (hashes.length != evidenceIdhashSet.size()) {
					log.debug("Hash collision when rehydrating evidence set for assembly " + getEvidenceID());
				}
				evidenceSetSize = hashes.length;
			}
		}
	}
	private static void setEvidenceSet(SAMRecord r, Collection<DirectedEvidence> evidence) {
		if (evidence != null && evidence.size() > 0) {
			int[] hashes = new int[evidence.size()];
			int i = 0;
			for (DirectedEvidence e : evidence) {
				hashes[i++] = evidenceIdHash.hashUnencodedChars(e.getEvidenceID()).asInt();
			}
			r.setAttribute(SamTags.ASSEMLBY_EVIDENCE_HASH, hashes);
		}
	}
	/**
	 * Determines whether the given record is part of the given assembly
	 *
	 * This method is a probabilistic method and it is possible for the record to return true
	 * when the record does not form part of the assembly breakend
	 *
	 * @param evidence
	 * @return true if the record is likely part of the breakend, false if definitely not
	 */
	public boolean isPartOfAssemblyBreakend(DirectedEvidence e) {
		ensureEvidenceSet();
		if (evidenceIdhashSet == null) return false;
		if (evidenceIdhashSet.contains(evidenceIdHash.hashUnencodedChars(e.getEvidenceID()).asInt())) {
			evidence.add(e);
			return true;
		} else {
			return false;
		}
	}
	public Collection<DirectedEvidence> getEvidence() {
		if (evidence == null) {
			return ImmutableList.of();
		}
		ensureEvidenceSet();
		if (evidenceIdhashSet != null) {
			if (evidence.size() < evidenceSetSize) {
				throw new IllegalStateException("Evidence not yet fully rehydrated");
			}
			if (evidence.size() > evidenceSetSize) {
				// Don't throw exception as the user can't actually do anything about this
				// Just continue with as correct results as we can manage
				log.warn("Hash collision has resulted in evidence being incorrectly attributed to assembly " + getEvidenceID());
			}
		}
		return evidence;
	}
	private void fixReadPair() {
		if (realignment.getReadUnmappedFlag()) {
			// SAMv1 S2.4
			realignment.setReferenceIndex(record.getReferenceIndex());
			realignment.setAlignmentStart(record.getAlignmentStart());
		}
		SAMRecordUtil.pairReads(this.record, this.realignment);
	}
	public SAMRecordAssemblyEvidence(
			Collection<DirectedEvidence> evidence,
			SAMFileHeader samFileHeader, BreakendSummary breakend,
			AssemblyEvidenceSource source,
			int anchoredBaseCount,
			byte[] baseCalls,
			byte[] baseQuals,
			Map<VcfAttributes, int[]> intListAttributes
			) {
		this(source, createAssemblySAMRecord(evidence, samFileHeader, source, breakend, anchoredBaseCount, baseCalls, baseQuals, intListAttributes), null);
		this.evidence = evidence;
	}
	private SAMRecord getPlaceholderRealignment() {
		SAMRecord placeholder = new SAMRecord(record.getHeader());
		placeholder.setReadUnmappedFlag(true);
		placeholder.setReadBases(getBreakendSequence());
		placeholder.setBaseQualities(getBreakendQuality());
		placeholder.setReadNegativeStrandFlag(false);
		return placeholder;
	}
	private static BreakendSummary calculateBreakendFromAlignmentCigar(Integer referenceIndex, int alignmentStart, Cigar cigar) {
		if (referenceIndex == null) return null;
		List<CigarElement> ce = cigar.getCigarElements();
		BreakendDirection direction = ce.get(0).getOperator() == CigarOperator.SOFT_CLIP ? BreakendDirection.Backward : BreakendDirection.Forward;
		int start = alignmentStart;
		int end = alignmentStart;
		int nonSoftClipCigarLength = 0;
		for (CigarElement e : ce) {
			if (e.getOperator() != CigarOperator.SOFT_CLIP) {
				nonSoftClipCigarLength += e.getLength();
			}
		}
		switch (ce.get(direction == BreakendDirection.Forward ? 0 : ce.size() - 1).getOperator()) {
			case M:
				assert(cigar.numCigarElements() == 2);
				if (direction == BreakendDirection.Forward) {
					assert(cigar.getCigarElement(0).getOperator() == CigarOperator.MATCH_OR_MISMATCH);
					assert(cigar.getCigarElement(1).getOperator() == CigarOperator.SOFT_CLIP);
					// breakend is at end
					int anchorLength = nonSoftClipCigarLength;
					start += anchorLength - 1;
					end += anchorLength - 1;
				} else {
					assert(cigar.numCigarElements() == 2);
					assert(cigar.getCigarElement(0).getOperator() == CigarOperator.SOFT_CLIP);
					assert(cigar.getCigarElement(1).getOperator() == CigarOperator.MATCH_OR_MISMATCH);
				}
				break;
			case X:
				if (direction == BreakendDirection.Forward) {
					end += nonSoftClipCigarLength - 1;
				} else {
					start -= nonSoftClipCigarLength - 1;
				}
				break;
			default:
				throw new IllegalArgumentException(cigar + "is not a valid assembly CIGAR");
		}
		return new BreakendSummary(referenceIndex, direction, start, end);
	}
	private boolean calculateIsBreakendExactFromCigar(Cigar cigar) {
		for (CigarElement ce : cigar.getCigarElements()) {
			if (ce.getOperator() == CigarOperator.M) return true;
			if (ce.getOperator() == CigarOperator.X) return false;
		}
		throw new IllegalArgumentException("Cannot determine breakend type from cigar " + cigar.toString());
	}
	private static SAMRecord createAssemblySAMRecord(
			Collection<DirectedEvidence> evidence,
			SAMFileHeader samFileHeader, AssemblyEvidenceSource source,
			BreakendSummary breakend,
			int anchoredBaseCount, byte[] baseCalls, byte[] baseQuals,
			Map<VcfAttributes, int[]> intListAttributes) {
		if (breakend instanceof BreakpointSummary) throw new IllegalArgumentException("Breakpoints not supported by this constructor");
		SAMRecord record = new SAMRecord(samFileHeader);
		record.setReferenceIndex(breakend.referenceIndex);
		record.setReadName(String.format("gridss_%s_%d_%d%s_%d",
				source.processContext.getDictionary().getSequence(breakend.referenceIndex).getSequenceName(),
				breakend.start, breakend.end,
				breakend.direction == BreakendDirection.Forward ? "f" : "b",
				Math.abs(Arrays.hashCode(baseCalls))) + anchoredBaseCount);
		setIntListAttributes(record, intListAttributes);
		if (anchoredBaseCount > 0) {
			record.setReadBases(baseCalls);
			record.setBaseQualities(baseQuals);
			if (breakend.start != breakend.end) {
				throw new IllegalArgumentException("Imprecisely anchored breakends not supported by this constructor");
			}
			if (breakend.direction == BreakendDirection.Forward) {
				record.setAlignmentStart(breakend.start - anchoredBaseCount + 1);
				record.setCigar(new Cigar(ImmutableList.of(
						new CigarElement(anchoredBaseCount, CigarOperator.MATCH_OR_MISMATCH),
						new CigarElement(baseCalls.length - anchoredBaseCount, CigarOperator.SOFT_CLIP))));
			} else {
				record.setAlignmentStart(breakend.start);
				record.setCigar(new Cigar(ImmutableList.of(
						new CigarElement(baseCalls.length - anchoredBaseCount, CigarOperator.SOFT_CLIP),
						new CigarElement(anchoredBaseCount, CigarOperator.MATCH_OR_MISMATCH))));
			}
		} else {
			// SAM spec requires at least one mapped base
			// to conform to this, we add a placeholder mismatched base to our read
			// in the furthest anchor position
			// and represent the breakend confidence interval as an N
			// interval anchored by Xs
			LinkedList<CigarElement> ce = new LinkedList<CigarElement>();
			int len = breakend.end - breakend.start + 1;
			int padBases;
			if (len <= 2) {
				ce.add(new CigarElement(len, CigarOperator.X));
				padBases = len;
			} else {
				ce.add(new CigarElement(1, CigarOperator.X));
				ce.add(new CigarElement(len - 2, CigarOperator.N));
				ce.add(new CigarElement(1, CigarOperator.X));
				padBases = 2;
			}
			if (breakend.direction == BreakendDirection.Forward) {
				record.setAlignmentStart(breakend.start);
				ce.addLast(new CigarElement(baseCalls.length, CigarOperator.SOFT_CLIP));
				record.setCigar(new Cigar(ce));
				record.setReadBases(Bytes.concat(PAD_BASES[padBases], baseCalls));
				record.setBaseQualities(Bytes.concat(PAD_QUALS[padBases], baseQuals));
			} else {
				record.setAlignmentStart(breakend.end);
				ce.addFirst(new CigarElement(baseCalls.length, CigarOperator.SOFT_CLIP));
				record.setCigar(new Cigar(ce));
				record.setReadBases(Bytes.concat(baseCalls, PAD_BASES[padBases]));
				record.setBaseQualities(Bytes.concat(baseQuals, PAD_QUALS[padBases]));
			}
		}
		setEvidenceSet(record, evidence);
		return record;
	}
	private static final byte[][] PAD_BASES = new byte[][] { new byte[] {}, new byte[] { 'N' }, new byte[] { 'N', 'N' } };
	private static final byte[][] PAD_QUALS = new byte[][] { new byte[] {}, new byte[] { 0 }, new byte[] { 0, 0 } };
	private static void setIntListAttributes(SAMRecord record, Map<VcfAttributes, int[]> intListAttributes) {
		if (intListAttributes == null) throw new IllegalArgumentException("Attribute list not supplied");
		// Stored in our MAPQ flag
		record.setMappingQuality(AttributeConverter.asInt(intListAttributes.get(VcfAttributes.ASSEMBLY_MAPQ_LOCAL_MAX), 0));
		intListAttributes.remove(VcfAttributes.ASSEMBLY_MAPQ_LOCAL_MAX);
		for (Entry<VcfAttributes, int[]> entry : intListAttributes.entrySet()) {
			VcfAttributes attr = entry.getKey();
			int[] value = entry.getValue();
			String samtag = attr.samTag();
			if (StringUtils.isBlank(samtag)) {
				throw new IllegalArgumentException(String.format("Attribute %s does not have a sam tag assigned", attr));
			}
			if (attr.infoHeader() != null && attr.infoHeader().isFixedCount() && attr.infoHeader().getCount() != value.length) {
				throw new IllegalArgumentException(String.format("Attribute %s has %d values - expected %d", attr, value.length, attr.infoHeader().getCount()));
			}
			if (value.length == 1) {
				record.setAttribute(attr.samTag(), value[0]);
			} else {
				record.setAttribute(attr.samTag(), value);
			}
		}
	}
	private int getAnchorLength() {
		CigarElement ce;
		if (breakend.direction == BreakendDirection.Forward) {
			ce = record.getCigar().getCigarElement(0);
		} else {
			ce = record.getCigar().getCigarElement(record.getCigarLength() - 1);
		}
		if (ce.getOperator() == CigarOperator.X) {
			return 0;
		}
		return ce.getLength();
	}
	private int getBreakendLength() {
		CigarElement ce;
		if (breakend.direction == BreakendDirection.Forward) {
			ce = record.getCigar().getCigarElement(record.getCigarLength() - 1);
			
		} else {
			ce = record.getCigar().getCigarElement(0);
		}
		assert(ce.getOperator() == CigarOperator.SOFT_CLIP);
		return ce.getLength();
	}
	private byte[] getAnchorBytes(byte[] data) {
		int len = getAnchorLength();
		if (breakend.direction == BreakendDirection.Forward) {
			return Arrays.copyOfRange(data, 0, len);
		} else {
			return Arrays.copyOfRange(data, data.length - len, data.length);
		}
	}
	private byte[] getBreakendBytes(byte[] data) {
		int len = getBreakendLength();
		if (breakend.direction == BreakendDirection.Forward) {
			return Arrays.copyOfRange(data, data.length - len, data.length);
		} else {
			return Arrays.copyOfRange(data, 0, len);
		}
	}
	@Override
	public BreakendSummary getBreakendSummary() {
		return breakend;
	}
	@Override
	public byte[] getBreakendSequence() {
		return getBreakendBytes(record.getReadBases());
	}
	@Override
	public byte[] getBreakendQuality() {
		return getBreakendBytes(record.getBaseQualities());
	}
	@Override
	public String getEvidenceID() {
		return record.getReadName();
	}
	@Override
	public AssemblyEvidenceSource getEvidenceSource() {
		return source;
	}
	@Override
	public int getLocalMapq() {
		return record.getMappingQuality();
	}
	@Override
	public int getLocalBaseLength() {
		return getAnchorLength();
	}
	@Override
	public int getLocalMaxBaseQual() {
		if (getAnchorLength() == 0) return 0;
		return UnsignedBytes.toInt(UnsignedBytes.max(getAssemblyAnchorQuals()));
	}
	@Override
	public int getLocalTotalBaseQual() {
		int sum = 0;
		byte[] data = getAssemblyAnchorQuals();
		for (int i = 0; i < data.length; i++) {
			sum += UnsignedBytes.toInt(data[i]);
		}
		return sum;
	}
	@Override
	public byte[] getAssemblySequence() {
		if (breakend.direction == BreakendDirection.Forward) {
			return Bytes.concat(getAssemblyAnchorSequence(), getBreakendSequence());
		} else {
			return Bytes.concat(getBreakendSequence(), getAssemblyAnchorSequence());
		}
	}
	@Override
	public byte[] getAssemblyAnchorSequence() {
		return getAnchorBytes(record.getReadBases());
	}
	public byte[] getAssemblyAnchorQuals() {
		return getAnchorBytes(record.getBaseQualities());
	}
	@Override
	public int getAssemblyAnchorLength() {
		return getAnchorLength();
	}
	@Override
	public int getAssemblyBaseCount(EvidenceSubset subset) {
		return AttributeConverter.asIntSumTN(record.getAttribute(VcfAttributes.ASSEMBLY_BASE_COUNT.samTag()), subset);
	}
	@Override
	public int getAssemblySupportCountReadPair(EvidenceSubset subset) {
		return AttributeConverter.asIntSumTN(record.getAttribute(VcfAttributes.ASSEMBLY_READPAIR_COUNT.samTag()), subset);
	}
	@Override
	public int getAssemblyReadPairLengthMax(EvidenceSubset subset) {
		return AttributeConverter.asIntMaxTN(record.getAttribute(VcfAttributes.ASSEMBLY_READPAIR_LENGTH_MAX.samTag()), subset);
	}
	@Override
	public int getAssemblySupportCountSoftClip(EvidenceSubset subset) {
		return AttributeConverter.asIntSumTN(record.getAttribute(VcfAttributes.ASSEMBLY_SOFTCLIP_COUNT.samTag()), subset);
	}
	@Override
	public int getAssemblySoftClipLengthTotal(EvidenceSubset subset) {
		return AttributeConverter.asIntSumTN(record.getAttribute(VcfAttributes.ASSEMBLY_SOFTCLIP_CLIPLENGTH_TOTAL.samTag()), subset);
	}
	@Override
	public int getAssemblySoftClipLengthMax(EvidenceSubset subset) {
		return AttributeConverter.asIntMaxTN(record.getAttribute(VcfAttributes.ASSEMBLY_SOFTCLIP_CLIPLENGTH_MAX.samTag()), subset);
	}
	@Override
	public boolean isAssemblyFiltered() {
		return StringUtils.isNotBlank((String)record.getAttribute(SamTags.ASSEMLBY_FILTERS));
	}
	@Override
	public void filterAssembly(VcfFilter reason) {
		String tag = (String)record.getAttribute(SamTags.ASSEMLBY_FILTERS);
		if (StringUtils.isBlank(tag)) {
			tag = reason.filter();
		} else if (!tag.contains(reason.filter())) {
			tag = tag + "," + reason.filter();
		}
		record.setAttribute(SamTags.ASSEMLBY_FILTERS, tag);
	}
	@Override
	public List<VcfFilter> getFilters() {
		List<VcfFilter> list = Lists.newArrayList();
		String filters = (String)record.getAttribute(SamTags.ASSEMLBY_FILTERS);
		if (!StringUtils.isEmpty(filters)) {
			for (String s : filters.split(",")) {
				list.add(VcfFilter.get(s));
			}
		}
		return list;
	}
	public SAMRecord getSAMRecord() {
		return record;
	}
	public SAMRecord getRemoteSAMRecord() {
		return realignment;
	}
	@Override
	public String toString() {
		return String.format("A  %s N=%s", getBreakendSummary(), getEvidenceID());
	}
	static final Ordering<SAMRecordAssemblyEvidence> BySAMCoordinate = new Ordering<SAMRecordAssemblyEvidence>() {
		private final SAMRecordCoordinateComparator cmp = new SAMRecordCoordinateComparator();
		@Override
		public int compare(SAMRecordAssemblyEvidence arg0, SAMRecordAssemblyEvidence arg1) {
			return cmp.compare(arg0.getSAMRecord(), arg1.getSAMRecord());
		}
	};
	@Override
	public boolean isBreakendExact() {
		return isExact;
	}
}
