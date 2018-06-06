package au.edu.wehi.idsv.sam;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.apache.commons.lang3.StringUtils;

import com.google.common.collect.Ordering;
import com.google.common.primitives.Ints;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.BreakendSummary;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.TextCigarCodec;

/**
 * SA Z Other canonical alignments in a chimeric alignment, formatted as a semicolon-delimited
 * list: (rname,pos,strand,CIGAR,mapQ,NM ;)+. Each element in the list represents a
 * part of the chimeric alignment. Conventionally, at a supplementary line, the first element
 * points to the primary line.
 * @author Daniel Cameron
 *
 */
public class ChimericAlignment {
	public final String rname;
	public final int pos;
	public final boolean isNegativeStrand;
	public final Cigar cigar;
	public final int mapq;
	public final Integer nm;
	public ChimericAlignment(String rname, int pos, boolean isNegativeStrand, Cigar cigar, int mapq, int nm) {
		this.rname = rname;
		this.pos = pos;
		this.isNegativeStrand = isNegativeStrand;
		this.cigar = cigar;
		this.mapq = mapq;
		this.nm = nm;
	}
	public ChimericAlignment(SAMRecord r) {
		this.rname = r.getReferenceName();
		this.pos = r.getAlignmentStart();
		this.isNegativeStrand = r.getReadNegativeStrandFlag();
		this.cigar = r.getCigar();
		this.mapq = r.getMappingQuality();
		this.nm = r.getIntegerAttribute(SAMTag.NM.name());
	}
	public ChimericAlignment(String str) {
		str = str.replace(';', ' ').trim();
		String[] splits = str.split(",");
		this.rname = splits[0];
		this.pos = Integer.parseInt(splits[1]);
		this.isNegativeStrand = "-".equals(splits[2]);
		this.cigar = TextCigarCodec.decode(splits[3]);
		this.mapq = Integer.parseInt(splits[4]);
		Integer nmParsed = null;
		try {
			if (splits.length >= 6) {
				nmParsed = Integer.parseInt(splits[5]);
			}
		} catch (NumberFormatException nfe) {
			// swallow and fall back to null
		}
		this.nm = nmParsed;
	}
	public static List<ChimericAlignment> getChimericAlignments(String sa) {
		if (StringUtils.isEmpty(sa)) return Collections.emptyList();
		List<ChimericAlignment> list = new ArrayList<ChimericAlignment>();
		String[] splits = sa.split(";");
		for (String s : splits) {
			if (!StringUtils.isEmpty(sa)) {
				list.add(new ChimericAlignment(s));
			}
		}
		return list;
	}
	public static List<ChimericAlignment> getChimericAlignments(SAMRecord r) {
		return getChimericAlignments(r.getStringAttribute(SAMTag.SA.name()));
	}
	private BreakendSummary startBreakend(SAMSequenceDictionary dict) {
		return new BreakendSummary(rnameToReferenceIndex(dict, rname), BreakendDirection.Backward, pos);
	}
	private BreakendSummary endBreakend(SAMSequenceDictionary dict) {
		int endpos = pos + cigar.getReferenceLength() - 1;
		return new BreakendSummary(rnameToReferenceIndex(dict, rname), BreakendDirection.Forward, endpos);
	}
	private static int rnameToReferenceIndex(SAMSequenceDictionary dict, String rname) {
		int index = dict.getSequenceIndex(rname);
		if (index < 0) {
			throw new IllegalArgumentException(String.format("Reference sequence %s not found in sequence dictionary.", rname));
		}
		return index;
	}
	public BreakendSummary successorBreakend(SAMSequenceDictionary dict) {
		return isNegativeStrand ? startBreakend(dict) : endBreakend(dict);
	}
	public BreakendSummary predecessorBreakend(SAMSequenceDictionary dict) {
		return isNegativeStrand ? endBreakend(dict) : startBreakend(dict);
	}
	/**
	 * Gets the read offset of the first aligned base.
	 * @return zero based offset from start of read based on fastq sequencing base order
	 */
	public int getFirstAlignedBaseReadOffset() {
		return isNegativeStrand ? CigarUtil.getEndClipLength(cigar.getCigarElements()) : CigarUtil.getStartClipLength(cigar.getCigarElements());
	}
	public int getLastAlignedBaseReadOffset() {
		return cigar.getReadLength() - 1 - (isNegativeStrand ? CigarUtil.getStartClipLength(cigar.getCigarElements()) : CigarUtil.getEndClipLength(cigar.getCigarElements()));
		
	}
	@Override
	public String toString() {
		return String.format("%s,%d,%s,%s,%d,%s", rname, pos, isNegativeStrand ? "-" : "+", cigar, mapq, nm == null ? "" : nm.toString());
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((cigar == null) ? 0 : cigar.hashCode());
		result = prime * result + (isNegativeStrand ? 1231 : 1237);
		result = prime * result + mapq;
		//result = prime * result + ((nm == null) ? 0 : nm.hashCode());
		result = prime * result + pos;
		result = prime * result + ((rname == null) ? 0 : rname.hashCode());
		return result;
	}
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		ChimericAlignment other = (ChimericAlignment) obj;
		if (cigar == null) {
			if (other.cigar != null)
				return false;
		} else if (!cigar.equals(other.cigar))
			return false;
		if (isNegativeStrand != other.isNegativeStrand)
			return false;
		if (mapq != other.mapq)
			return false;
		//if (nm == null) {
		//	if (other.nm != null)
		//		return false;
		//} else if (!nm.equals(other.nm))
		//	return false;
		if (pos != other.pos)
			return false;
		if (rname == null) {
			if (other.rname != null)
				return false;
		} else if (!rname.equals(other.rname))
			return false;
		return true;
	}
	public static final Ordering<ChimericAlignment> ByReadOffset = new Ordering<ChimericAlignment>() {
		@Override
		public int compare(ChimericAlignment left, ChimericAlignment right) {
			return Ints.compare(left.getFirstAlignedBaseReadOffset(), right.getFirstAlignedBaseReadOffset());
		}
	};
}
