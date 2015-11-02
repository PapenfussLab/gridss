package au.edu.wehi.idsv.sam;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.lang3.StringUtils;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;
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
	public final int nm;
	public ChimericAlignment(String rname, int pos, boolean isNegativeStrand, Cigar cigar, int mapq, int nm) {
		this.rname = rname;
		this.pos = pos;
		this.isNegativeStrand = isNegativeStrand;
		this.cigar = cigar;
		this.mapq = mapq;
		this.nm = nm;
	}
	public ChimericAlignment(String str) {
		str = str.replace(';', ' ').trim();
		String[] splits = str.split(",");
		this.rname = splits[0];
		this.pos = Integer.parseInt(splits[1]);
		this.isNegativeStrand = "-".equals(splits[2]);
		this.cigar = TextCigarCodec.decode(splits[3]);
		this.mapq =Integer.parseInt(splits[4]);
		this.nm = Integer.parseInt(splits[5]);
	}
	public static List<ChimericAlignment> getChimericAlignments(String sa) {
		List<ChimericAlignment> list = new ArrayList<ChimericAlignment>();
		if (sa != null) {
			String[] splits = sa.split(";");
			for (String s : splits) {
				if (!StringUtils.isEmpty(sa)) {
					list.add(new ChimericAlignment(s));
				}
			}
		}
		return list;
	}
	public static List<ChimericAlignment> getChimericAlignments(SAMRecord r) {
		return getChimericAlignments(r.getStringAttribute("SA"));
	}
}
