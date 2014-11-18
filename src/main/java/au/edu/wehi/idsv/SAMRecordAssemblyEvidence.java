package au.edu.wehi.idsv;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.apache.commons.lang3.StringUtils;

import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.sam.SamTags;
import au.edu.wehi.idsv.vcf.VcfAttributes;
import au.edu.wehi.idsv.vcf.VcfFilter;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Ordering;
import com.google.common.primitives.Bytes;
import com.google.common.primitives.UnsignedBytes;

public class SAMRecordAssemblyEvidence implements AssemblyEvidence {
	private final SAMRecord record;
	private final SAMRecord realignment;
	private final AssemblyEvidenceSource source;
	private final BreakendSummary breakend;
	public SAMRecordAssemblyEvidence(AssemblyEvidenceSource source, SAMRecord assembly, SAMRecord realignment) {
		this.source = source;
		this.record = assembly;
		this.breakend = calculateBreakendFromAlignmentCigar(assembly.getReferenceIndex(), assembly.getAlignmentStart(), assembly.getCigar());
		this.realignment = realignment == null ? getPlaceholderRealignment() : realignment;
		fixReadPair();
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
			SAMFileHeader samFileHeader, BreakendSummary breakend,
			AssemblyEvidenceSource source,
			int anchoredBaseCount,
			byte[] baseCalls,
			byte[] baseQuals,
			Map<VcfAttributes, int[]> intListAttributes
			) {
		this(source, createAssemblySAMRecord(samFileHeader, source, breakend, anchoredBaseCount, baseCalls, baseQuals, intListAttributes), null);
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
		BreakendDirection direction = cigar.getCigarElement(0).getOperator() == CigarOperator.SOFT_CLIP ? BreakendDirection.Backward : BreakendDirection.Forward;
		int start = alignmentStart;
		int end = alignmentStart;
		if (cigar.numCigarElements() > 2) {
			int paddingNs = cigar.getCigarElement(1).getLength();
			if (direction == BreakendDirection.Forward) {
				end += paddingNs;
				assert(cigar.numCigarElements() == 3);
				assert(cigar.getCigarElement(0).getOperator() == CigarOperator.X);
				assert(cigar.getCigarElement(1).getOperator() == CigarOperator.SKIPPED_REGION);
				assert(cigar.getCigarElement(2).getOperator() == CigarOperator.SOFT_CLIP);
			} else {
				start -= paddingNs;
				assert(cigar.numCigarElements() == 3);
				assert(cigar.getCigarElement(0).getOperator() == CigarOperator.SOFT_CLIP);
				assert(cigar.getCigarElement(1).getOperator() == CigarOperator.SKIPPED_REGION);
				assert(cigar.getCigarElement(2).getOperator() == CigarOperator.X);
			}
		} else {
			if (direction == BreakendDirection.Forward) {
				assert(cigar.numCigarElements() == 2);
				assert(cigar.getCigarElement(0).getOperator() == CigarOperator.MATCH_OR_MISMATCH);
				assert(cigar.getCigarElement(1).getOperator() == CigarOperator.SOFT_CLIP);
			} else {
				assert(cigar.numCigarElements() == 2);
				assert(cigar.getCigarElement(0).getOperator() == CigarOperator.SOFT_CLIP);
				assert(cigar.getCigarElement(1).getOperator() == CigarOperator.MATCH_OR_MISMATCH);
			}
		}
		return new BreakendSummary(referenceIndex, direction, start, end);
	}
	private static SAMRecord createAssemblySAMRecord(
			SAMFileHeader samFileHeader, AssemblyEvidenceSource source,
			BreakendSummary breakend,
			int anchoredBaseCount, byte[] baseCalls, byte[] baseQuals,
			Map<VcfAttributes, int[]> intListAttributes) {
		if (breakend instanceof BreakpointSummary) throw new IllegalArgumentException("Breakpoints not supported by this constructor");
		SAMRecord record = new SAMRecord(samFileHeader);
		record.setReferenceIndex(breakend.referenceIndex);
		record.setReadName(String.format("assembly_%s_%d_%d_%d",
				source.processContext.getDictionary().getSequence(breakend.referenceIndex).getSequenceName(),
				breakend.start, breakend.end,
				Math.abs(Arrays.hashCode(baseCalls))));
		setIntListAttributes(record, intListAttributes);
		if (anchoredBaseCount > 0) {
			record.setReadBases(baseCalls);
			record.setBaseQualities(baseQuals);
			if (breakend.start != breakend.end) throw new IllegalArgumentException("Inexact anchored breakends not supported by this constructor");
			if (breakend.direction == BreakendDirection.Forward) {
				record.setAlignmentStart(breakend.start - anchoredBaseCount - 1);
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
			// and represent the breakend confidence interval with a placeholder
			// cigar N
			if (breakend.direction == BreakendDirection.Forward) {
				record.setAlignmentStart(breakend.start);
				record.setCigar(new Cigar(ImmutableList.of(
						new CigarElement(1, CigarOperator.X),
						new CigarElement(breakend.end - breakend.start, CigarOperator.SKIPPED_REGION),
						new CigarElement(baseCalls.length, CigarOperator.SOFT_CLIP))));
				record.setReadBases(Bytes.concat(new byte[] { 'N' }, baseCalls));
				record.setBaseQualities(Bytes.concat(new byte[] { 0 }, baseQuals));
			} else {
				record.setAlignmentStart(breakend.end);
				record.setCigar(new Cigar(ImmutableList.of(
						new CigarElement(baseCalls.length, CigarOperator.SOFT_CLIP),
						new CigarElement(breakend.end - breakend.start, CigarOperator.SKIPPED_REGION),
						new CigarElement(1, CigarOperator.X))));
				record.setReadBases(Bytes.concat(baseCalls, new byte[] { 'N' }));
				record.setBaseQualities(Bytes.concat(baseQuals, new byte[] { 0 }));
			}
		}
		return record;
	}
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
		List<VcfFilter> list = new ArrayList<>();
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
	static final Ordering<SAMRecordAssemblyEvidence> BySAMCoordinate = new Ordering<SAMRecordAssemblyEvidence>() {
		@Override
		public int compare(SAMRecordAssemblyEvidence arg0, SAMRecordAssemblyEvidence arg1) {
			return new SAMRecordCoordinateComparator().compare(arg0.getSAMRecord(), arg1.getSAMRecord());
		}
	};
}
