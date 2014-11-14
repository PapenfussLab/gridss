package au.edu.wehi.idsv;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.apache.commons.lang3.StringUtils;

import au.edu.wehi.idsv.sam.SamTags;
import au.edu.wehi.idsv.vcf.VcfAttributes;
import au.edu.wehi.idsv.vcf.VcfFilter;

import com.google.common.collect.ImmutableList;
import com.google.common.primitives.Bytes;
import com.google.common.primitives.UnsignedBytes;

public class SAMRecordAssemblyEvidence implements AssemblyEvidence {
	private SAMRecord record;
	private AssemblyEvidenceSource source;
	private BreakendSummary breakend;
	public SAMRecordAssemblyEvidence(
			BreakendSummary breakend,
			AssemblyEvidenceSource source,
			int anchoredBaseCount,
			byte[] baseCalls,
			byte[] baseQuals,
			Map<VcfAttributes, int[]> intListAttributes
			) {
		if (breakend instanceof BreakpointSummary) throw new IllegalArgumentException("Breakpoints not supported by this constructor");
		this.breakend = breakend;
		this.source = source;
		this.record = new SAMRecord(null);
		this.record.setReferenceIndex(breakend.referenceIndex);
		this.record.setReadName(String.format("assembly_%s_%d_%d_%d",
				source.processContext.getDictionary().getSequence(breakend.referenceIndex).getSequenceName(),
				breakend.start, breakend.end,
				Math.abs(Arrays.hashCode(baseCalls))));
		setIntListAttributes(record, intListAttributes);
		if (anchoredBaseCount > 0) {
			this.record.setReadBases(baseCalls);
			this.record.setBaseQualities(baseQuals);
			if (breakend.start != breakend.end) throw new IllegalArgumentException("Inexact anchored breakends not supported by this constructor");
			if (breakend.direction == BreakendDirection.Forward) {
				this.record.setAlignmentStart(breakend.start - anchoredBaseCount - 1);
				this.record.setCigar(new Cigar(ImmutableList.of(
						new CigarElement(anchoredBaseCount, CigarOperator.MATCH_OR_MISMATCH),
						new CigarElement(baseCalls.length - anchoredBaseCount, CigarOperator.SOFT_CLIP))));
			} else {
				this.record.setAlignmentStart(breakend.start);
				this.record.setCigar(new Cigar(ImmutableList.of(
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
				this.record.setAlignmentStart(breakend.start);
				this.record.setCigar(new Cigar(ImmutableList.of(
						new CigarElement(1, CigarOperator.X),
						new CigarElement(breakend.end - breakend.start - 1, CigarOperator.SKIPPED_REGION),
						new CigarElement(baseCalls.length, CigarOperator.SOFT_CLIP))));
				this.record.setReadBases(Bytes.concat(new byte[] { 'N' }, baseCalls));
				this.record.setBaseQualities(Bytes.concat(new byte[] { 0 }, baseQuals));
			}
			this.record.setAlignmentStart(breakend.end);
			this.record.setCigar(new Cigar(ImmutableList.of(
					new CigarElement(baseCalls.length, CigarOperator.SOFT_CLIP),
					new CigarElement(breakend.end - breakend.start - 1, CigarOperator.SKIPPED_REGION),
					new CigarElement(1, CigarOperator.X))));
			this.record.setReadBases(Bytes.concat(baseCalls, new byte[] { 'N' }));
			this.record.setBaseQualities(Bytes.concat(baseQuals, new byte[] { 0 }));
		}
	}
	private static void setIntListAttributes(SAMRecord record, Map<VcfAttributes, int[]> intListAttributes) {
		for (Entry<VcfAttributes, int[]> entry : intListAttributes.entrySet()) {
			VcfAttributes attr = entry.getKey();
			int[] value = entry.getValue();
			String samtag = attr.samTag();
			if (StringUtils.isBlank(samtag)) {
				throw new IllegalArgumentException(String.format("Attribute %s does not have a sam tag assigned", samtag));
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
		// Stored in our MAPQ flag
		record.setAttribute(VcfAttributes.ASSEMBLY_MAPQ_LOCAL_MAX.samTag(), AttributeConverter.asInt(intListAttributes.get(VcfAttributes.ASSEMBLY_MAPQ_LOCAL_MAX), 0));
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
	public EvidenceSource getEvidenceSource() {
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
	public int getLocalBaseCount() {
		throw new RuntimeException("NYI");
		//return getAnchorLength();
	}
	@Override
	public int getLocalMaxBaseQual() {
		if (getAnchorLength() != 0) return 0;
		return UnsignedBytes.toInt(UnsignedBytes.max(getAssemblyAnchorSequence()));
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
		if (!tag.contains(reason.filter())) {
			if (StringUtils.isBlank(tag)) {
				tag = reason.filter();
			} else {
				tag = tag + "," + reason.filter();
			}
			record.setAttribute(SamTags.ASSEMLBY_FILTERS, tag);
		}
	}
}
