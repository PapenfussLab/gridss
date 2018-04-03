package au.edu.wehi.idsv;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import org.apache.commons.lang3.ArrayUtils;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;

import au.edu.wehi.idsv.sam.ChimericAlignment;
import au.edu.wehi.idsv.sam.CigarUtil;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecord.SAMTagAndValue;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.SequenceUtil;

/**
 * Helper class for split read identification
 *
 */
public class SplitReadIdentificationHelper {
	private static final char SEPARATOR = '#';
	private static Comparator<SAMRecord> ByFirstAlignedBaseReadOffset = Comparator.comparing(r -> SAMRecordUtil.getFirstAlignedBaseReadOffset(r));
	/**
	 * Extract the unaligned portions of the read requiring realignment to identify split reads
	 * 
	 * @param r alignment record
	 * @param recordIsPartialAlignment true if the record is the result of aligning FastqRecords
	 * from a previous call to getSplitReadRealignments() 
	 * @return bases requiring alignment to identify split reads
	 */
	public static FastqRecord getFullRealignment(SAMRecord r, EvidenceIdentifierGenerator eidgen) {
		String name = eidgen.getAlignmentUniqueName(r);
		byte[] seq = r.getReadBases();
		byte[] qual = r.getBaseQualities();
		if (r.getReadNegativeStrandFlag()) {
			SequenceUtil.reverseComplement(seq);
			ArrayUtils.reverse(qual);
		}
		FastqRecord fastq = new FastqRecord(
				getSplitAlignmentFastqName(name, 0),
				new String(seq, StandardCharsets.US_ASCII),
				"",
				SAMUtils.phredToFastq(qual));
		return fastq;
	}
	public static List<FastqRecord> getSplitReadRealignments(SAMRecord r, boolean recordIsPartialAlignment, EvidenceIdentifierGenerator eidgen) {
		int startClipLength = SAMRecordUtil.getStartSoftClipLength(r);
		int endClipLength = SAMRecordUtil.getEndSoftClipLength(r);
		if (startClipLength + endClipLength == 0 || r.getReadUnmappedFlag()) {
			return Collections.emptyList();
		}
		int offset;
		String name;
		if (recordIsPartialAlignment) {
			offset = getRealignmentFirstAlignedBaseReadOffset(r);
			name = getOriginatingAlignmentUniqueName(r);
		} else {
			offset = 0;
			name = eidgen.getAlignmentUniqueName(r);
		}
		List<FastqRecord> list = new ArrayList<>(2);
		if (startClipLength > 0) {
			byte[] startbases = Arrays.copyOfRange(r.getReadBases(), 0, startClipLength);
			byte[] startqual  = Arrays.copyOfRange(r.getBaseQualities(), 0, startClipLength);
			int startOffset = offset;
			if (r.getReadNegativeStrandFlag()) {
				SequenceUtil.reverseComplement(startbases);
				ArrayUtils.reverse(startqual);
				// 5432
				// SS
				startOffset = startOffset + r.getReadLength() - startClipLength;  
			}
			FastqRecord start = new FastqRecord(
					getSplitAlignmentFastqName(name, startOffset),
					new String(startbases, StandardCharsets.US_ASCII),
					"",
					SAMUtils.phredToFastq(startqual));
			list.add(start);
		}
		if (endClipLength > 0) {
			byte[] endbases = Arrays.copyOfRange(r.getReadBases(), r.getReadLength() - endClipLength, r.getReadLength());
			byte[] endqual  = Arrays.copyOfRange(r.getBaseQualities(), r.getReadLength() - endClipLength, r.getReadLength());
			int endOffset = offset + r.getReadLength() - endClipLength;  
			if (r.getReadNegativeStrandFlag()) {
				SequenceUtil.reverseComplement(endbases);
				ArrayUtils.reverse(endqual);
				// 5432
				// SS
				endOffset = offset;  
			}
			FastqRecord end = new FastqRecord(
					getSplitAlignmentFastqName(name, endOffset),
					new String(endbases, StandardCharsets.US_ASCII),
					"",
					SAMUtils.phredToFastq(endqual));
			list.add(end);
		}
		return list;
    }
	/**
	 * Updates the given record and realignments to represent a split read alignment
	 * @param record initial read alignment
	 * @param salist realignments created from recursive alignment of the output of getSplitReadRealignments()
	 */
	public static void convertToSplitRead(SAMRecord record, List<SAMRecord> salist) {
		if (salist.size() == 0) return;
		List<SAMRecord> alignments = new ArrayList<>(1 + salist.size());
		for (SAMRecord r : salist) {
			if (!r.getReadUnmappedFlag()) {
				alignments.add(r);
			}
		}
		if (alignments.size() == 0) return;
		// ok, so we have an actual split read alignment
		unclip(record, alignments);
		for (SAMRecord r : alignments) {
			convertToSupplementaryAlignmentRecord(record, r);
		}
		writeSA(record, alignments);
	}
	/**
	 * Writes the SA SAM tag for split read alignments
	 * @param record primary alignment record
	 * @param alignments supplementary alignments
	 */
	private static void writeSA(SAMRecord record, List<SAMRecord> alignments) {
		alignments.sort(ByFirstAlignedBaseReadOffset);
		List<String> satags = alignments
				.stream()
				.map(r -> new ChimericAlignment(r).toString())
				.collect(Collectors.toList());
		record.setAttribute(SAMTag.SA.name(), satags.stream().collect(Collectors.joining(";")));
		// Primary alignment is first record as per SAMTags specs
		satags.add(0, new ChimericAlignment(record).toString());
		for (int i = 0; i < alignments.size(); i++) {
			final String currentsatag = satags.get(i + 1);
			String tag = satags.stream()
					.filter(s -> s != currentsatag) // can use reference equality since we're comparing to ourself
					.collect(Collectors.joining(";"));
			alignments.get(i).setAttribute(SAMTag.SA.name(), tag);
		}
	}
	/**
	 * Converts partial alignments to soft clipped alignments of the entire read
	 * @param record initial alignment containing the full record
	 * @param alignments read generated from recursive split read alignment using getSplitReadRealignments() 
	 */
	private static void unclip(SAMRecord record, List<SAMRecord> alignments) {
		byte[] readbases = record.getReadBases();
		byte[] quals = record.getBaseQualities();
		byte[] nreadbases = readbases.clone();
		byte[] nquals = quals.clone();
		SequenceUtil.reverseComplement(nreadbases);
		ArrayUtils.reverse(nquals);
		if (record.getReadNegativeStrandFlag()) {
			byte[] tmp = readbases;
			readbases = nreadbases;
			nreadbases = tmp;
			tmp = quals;
			quals = nquals;
			nquals = tmp;
		}
		// convert to chimeric alignments
		for (SAMRecord r : alignments) {
			int offset = getRealignmentFirstAlignedBaseReadOffset(r);
			int prepad = offset;
			int postpad = record.getReadLength() - offset - r.getReadLength();
			if (r.getReadNegativeStrandFlag()) {
				r.setReadBases(nreadbases);
				r.setBaseQualities(nquals);
			} else {
				r.setReadBases(readbases);
				r.setBaseQualities(quals);
			}
			padCigar(r, prepad, postpad);
		}
	}
	private static void padCigar(SAMRecord record, int pre, int post) {
		List<CigarElement> cigar = Lists.newArrayList(record.getCigar().getCigarElements());
		if (record.getReadNegativeStrandFlag()) {
			cigar.add(0, new CigarElement(post, CigarOperator.SOFT_CLIP));
			cigar.add(new CigarElement(pre, CigarOperator.SOFT_CLIP));
		} else {
			cigar.add(0, new CigarElement(pre, CigarOperator.SOFT_CLIP));
			cigar.add(new CigarElement(post, CigarOperator.SOFT_CLIP));
		}
		cigar = CigarUtil.clean(cigar, false);
		record.setCigar(new Cigar(cigar));
	}
	private static void convertToSupplementaryAlignmentRecord(SAMRecord primary, SAMRecord supplementary) {
		supplementary.setReadName(primary.getReadName());
		supplementary.setMateAlignmentStart(primary.getMateAlignmentStart());
		supplementary.setMateReferenceIndex(primary.getMateReferenceIndex());
		supplementary.setReadPairedFlag(primary.getReadPairedFlag());
		if (primary.getReadPairedFlag()) {
			supplementary.setProperPairFlag(primary.getProperPairFlag());
			supplementary.setFirstOfPairFlag(primary.getFirstOfPairFlag());
			supplementary.setSecondOfPairFlag(primary.getSecondOfPairFlag());
			supplementary.setMateUnmappedFlag(primary.getMateUnmappedFlag());
			supplementary.setMateNegativeStrandFlag(primary.getMateNegativeStrandFlag());
		}
		supplementary.setSecondaryAlignment(primary.isSecondaryAlignment());
		supplementary.setDuplicateReadFlag(primary.getDuplicateReadFlag());
		supplementary.setReadFailsVendorQualityCheckFlag(primary.getReadFailsVendorQualityCheckFlag());
		supplementary.setSupplementaryAlignmentFlag(true);
		// attributes not already set by the supplementary alignment 
		for (SAMTagAndValue attr : primary.getAttributes()) {
			if (supplementary.getAttribute(attr.tag) == null) {
				supplementary.setAttribute(attr.tag, attr.value);
			}
		}
	}
	public static String getOriginatingAlignmentUniqueName(SAMRecord splitread) {
		String name = splitread.getReadName();
		name = name.substring(0, name.lastIndexOf(SEPARATOR));
		return name;
	}
	public static int getRealignmentFirstAlignedBaseReadOffset(SAMRecord splitread) {
		String name = splitread.getReadName();
		name = name.substring(name.lastIndexOf(SEPARATOR) + 1);
		return Integer.parseInt(name);
	}
	public static String getSplitAlignmentFastqName(String alignmentUniqueName, int firstAlignedBaseReadOffset) {
		String fastqid = alignmentUniqueName + SEPARATOR + Integer.toString(firstAlignedBaseReadOffset);
		return fastqid;
	}
	/**
	 * Replaces the primary alignment with one of the given alignments.
	 * Replacement preference is given to alignments that partially overlap the originatingRecord alignment
	 * @param originatingRecord
	 * @param realignments
	 * @return alignment record that the originatingRecord alignment was replaced with.
	 */
	public static SAMRecord replaceAlignment(SAMRecord originatingRecord, List<SAMRecord> realignments) {
		// Find alignment that best matches the originatingRecord
		SAMRecord newPrimary = null;
		int maxOverlap = 0; 
		for (SAMRecord r : realignments) {
			int overlap = SAMRecordUtil.overlappingBases(originatingRecord, r);
			if (overlap > maxOverlap) {
				maxOverlap = overlap;
				newPrimary = r;
			}
		}
		// otherwise just get the first alignment returned by the aligner 
		if (newPrimary == null && realignments.size() > 0) {
			newPrimary = realignments.get(0);
		}
		replaceAlignment(originatingRecord, newPrimary);
		return newPrimary;
	}
	private static void replaceAlignment(SAMRecord originatingRecord, SAMRecord newPrimary) {
		originatingRecord.setAttribute("OA", new ChimericAlignment(originatingRecord).toString());
		if (newPrimary == null || newPrimary.getReadUnmappedFlag()) {
			originatingRecord.setReadUnmappedFlag(true);
		} else {
			unclip(originatingRecord, ImmutableList.of(newPrimary));
			originatingRecord.setReadUnmappedFlag(newPrimary.getReadUnmappedFlag());
			originatingRecord.setAlignmentStart(newPrimary.getAlignmentStart());
			originatingRecord.setReferenceIndex(newPrimary.getReferenceIndex());
			originatingRecord.setCigar(newPrimary.getCigar());
			if (originatingRecord.getReadNegativeStrandFlag() != newPrimary.getReadNegativeStrandFlag()) {
				byte[] b = originatingRecord.getReadBases();
				SequenceUtil.reverseComplement(b);
				originatingRecord.setReadBases(b);
				byte[] q = originatingRecord.getBaseQualities();
				ArrayUtils.reverse(q);
				originatingRecord.setBaseQualities(q);
			}
			originatingRecord.setReadNegativeStrandFlag(newPrimary.getReadNegativeStrandFlag());
		}
	}
}
