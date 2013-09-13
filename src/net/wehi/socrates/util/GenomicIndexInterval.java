/**
 * 
 */
package net.wehi.socrates.util;

//import net.sf.samtools.SAMRecord;

/**
 * 
 */

import java.io.Serializable;

import org.apache.commons.lang3.builder.HashCodeBuilder;
/**
 * @author Arthur Hsu
 * Created: Sep 24, 2012
 *
 */
public class GenomicIndexInterval implements Serializable, Comparable<GenomicIndexInterval> {
	private static final long serialVersionUID = 61393454553L; 
	public int chromIdx;
	public int start, end;
	public byte strand;
	
	/**
	 * Clone constructor.
	 * @param gi
	 */
	public GenomicIndexInterval(GenomicIndexInterval gi) {
		chromIdx = gi.chromIdx;
		start = gi.start;
		end = gi.end;
		strand = gi.strand;
	}
	
	/**
	 * Constructor. Defaults to unstranded interval.
	 * @param chrom
	 * @param start
	 * @param end
	 */
	public GenomicIndexInterval(int chromIdx, int start, int end) {
		this(chromIdx, start, end, '.');
	}

	/**
	 * Constructor. Specifying all parameters.
	 * @param chrom
	 * @param start
	 * @param end
	 * @param strand
	 */
	public GenomicIndexInterval(int chromIdx, int start, int end, char strand) {
		this.chromIdx = chromIdx;
		int s = (start<0 || end<0) ? Math.max(start, end) : Math.min(start,end);
		int e = (start<0 || end<0) ? s : Math.max(start, end);
		this.start = s;
		this.end = e;
		
		if (strand==' ') { // default is to have strand determined by coordinate
			this.strand = (s!=start) ? (byte)-1 : (byte)1;
		} else {
			this.strand = (strand=='+') ? (byte)1 : ((strand=='-') ? (byte)-1 : 0);
		}
	}

	/**
 	 * Provide ordering of two GenomicInterval classes. First order by chromosome (lexigraphically), then by coordinate. 
	 * @param other
	 * @return -1 if less than; 0 if equal; +1 if greater
	 */
	public int compareTo(GenomicIndexInterval other) {
		if (this.chromIdx != other.chromIdx) return this.chromIdx-other.chromIdx;
		if (this.start != other.start) return this.start-other.start;
    	if (this.strand != other.strand) return other.strand-this.strand;
    	if (this.end != other.end) return this.end-other.end;
    	return 0;
	}
	
	@Override
	public boolean equals(Object other) {
		GenomicIndexInterval o = (GenomicIndexInterval)other;
		if (o.chromIdx==chromIdx && o.start==start && o.end==end && o.strand==strand) return true;
		return false;
	}
	
	@Override
	public int hashCode() { return new HashCodeBuilder(13, 17).append(chromIdx).append(start).append(end).append(strand).toHashCode(); }
	
	public boolean adjoins(GenomicIndexInterval other) {
		if (chromIdx!=other.chromIdx) return false;
		if (intersects(other)) return false;
		if (Math.abs(other.start-end)==1 || Math.abs(start-other.end)==1) return true;
		else return false;
	}
	
	/***
	 * Tests if the genomic interval intersects another.
	 * @param other
	 * @return true if intervals intersect
	 */
	public boolean intersects(GenomicIndexInterval other) {
		return this.intersects(other, 0);
	}
	
	/**
	 * Tests if the genomic interval intersects another with a flanking region. 
	 * @param other
	 * @param flank
	 * @return true if intervals intersect
	 */
	public boolean intersects(GenomicIndexInterval other, int flank) {
		if (this.chromIdx != other.chromIdx) return false;
		if ( (this.strand != 0 && other.strand != 0) && this.strand != other.strand) return false;
		if (other.start <= this.end+flank && other.end >= this.start-flank) return true;
		else return false;
	}
	
//	public boolean intersects(SAMRecord aln) {
//		if (this.chromIdx != aln.getReferenceIndex()) return false;
//		if (this.strand != 0) {
//			byte alnStrand = aln.getReadNegativeStrandFlag() ? (byte)-1 : (byte)1;
//			if ( alnStrand != this.strand ) return false;
//		}
//		if (aln.getAlignmentStart() <= this.end && aln.getAlignmentEnd() >= this.start) return true;
//		else return false;
//	}

	/***
	 * Tests if the genomic interval intersects another, ignoring strand.
	 * @param other
	 * @return true if intervals intersect regardless of strand
	 */
	public boolean intersectsIgnoreStrand(GenomicIndexInterval other) {
		return this.intersectsIgnoreStrand(other, 0);
	}
	
	/**
	 * Tests if the genomic interval intersects another with a flanking region - Ignore strand. 
	 * @param other
	 * @param flank
	 * @return true if intervals intersect regardless of strand
	 */
	public boolean intersectsIgnoreStrand(GenomicIndexInterval other, int flank) {
		if (this.chromIdx != other.chromIdx) return false;
		if (other.start <= this.end+flank && other.end >= this.start-flank) return true;
		else return false;
	}

	/***
	 * Tests if the genomic interval intersects another, having the opposite strand.
	 * @param other
	 * @return true if intervals of opposite strands intersect
	 */
	public boolean intersectsOppositeStrand(GenomicIndexInterval other) {
		return this.intersectsOppositeStrand(other, 0);
	}
	
	/**
	 * Tests if the genomic interval intersects another with a flanking region - having the opposite strand. 
	 * @param other
	 * @param flank
	 * @return true if intervals of opposite strands intersect
	 */
	public boolean intersectsOppositeStrand(GenomicIndexInterval other, int flank) {
		if (this.chromIdx != other.chromIdx) return false;
		if ( (this.strand != 0 && other.strand != 0) && this.strand == other.strand) return false;
		if (other.start <= this.end+flank && other.end >= this.start-flank) return true;
		else return false;
	}

	/**
	 * Tests if the genomic interval fully contains the other.
	 * @param other
	 * @return true if this interval fully contains the other
	 */
	public boolean contains(GenomicIndexInterval other) {
		if (this.chromIdx == other.chromIdx) return false;
		if ( (this.strand != 0 && other.strand != 0) && this.strand != other.strand) return false;
		return this.contains(other, 0);
	}

	/**
	 * Tests if the genomic interval fully contains the other with a flanking region.
	 * @param other
	 * @param flank
	 * @return true if this interval fully contains the other
	 */
	public boolean contains(GenomicIndexInterval other, int flank) {
		if (this.chromIdx != other.chromIdx) return false;
		if ( (this.strand != 0 && other.strand != 0) && this.strand != other.strand) return false;
		if (start-flank <= other.start && other.end <= end+flank) return true;
		return false;
	}
	
	public char getStrand() {
		if (strand==-1) return '-';
		else if (strand==1) return '+';
		else return '.';
	}
	
	/**
	 * Convert to String
	 * @param header
	 * @return Outputs the coordinate in the form of [chrom]:[start]-[end]
	 */
	public String toString(SAMFileInfo header) {
		switch (strand) {
		case -1:
			return header.getSequenceName(chromIdx) + ":" + start + "-" + end + "\t-";
		case 1:
			return header.getSequenceName(chromIdx) + ":" + start + "-" + end + "\t+";
		default:
			return header.getSequenceName(chromIdx) + ":" + start + "-" + end;
		}
	}
	
	
	/**
	 * Parse a string of genomic coordinate into usable form.
	 * @param region String in the format of, e.g., "chr1:1000000-1001000".
	 * @return Unstranded GenomicIndexInterval
	 */
	public static GenomicIndexInterval parse(String region, SAMFileInfo header) {
		// strip all "," in numeral coordinate
		region = region.replace(",", "");
		
		// split coordinate
		String[] tokens = region.split(":");
		
		String seq = null;
		int sid = -1;
		int start = -1, end = -1;
		if (tokens.length==1) {
			seq = new String(tokens[0]);
			if (header.getSequenceIndex(seq)==-1) {
				System.err.println("Sequence " + seq + " not found.");
				System.exit(1);
			}
			sid = header.getSequenceIndex(seq);
			start = 1;
			end = header.sequenceLengths.get(seq);
		} else if (tokens.length==2){
			seq = new String(tokens[0]);
			String[] posToken = tokens[1].split("-");
			if (posToken.length==2) {
				start = Math.max(1,Integer.parseInt(posToken[0]));
				end = Math.min(header.sequenceLengths.get(seq),Integer.parseInt(posToken[1]));
			} else if (posToken.length==1){
				start = Math.max(1,Integer.parseInt(posToken[0]));;
			} else {
				System.err.println("Invalid region format.");
				System.exit(1);
			}
			if (header.getSequenceIndex(seq)==-1) {
				System.err.println("Sequence " + seq + " not found.");
				System.exit(1);
			}
			sid = header.getSequenceIndex(seq);
		} else {
			System.err.println("Invalid region format.");
			System.exit(1);
		}
		return new GenomicIndexInterval( sid, start, end, '.' );
	}

	/**
	 * Create a new instance of GenomicInterval that is shifted from this one by "offset" (i.e. negative number to move upstream).
	 * @param offset
	 * @return shifted interval
	 */
	public GenomicIndexInterval createShiftedInterval(int offset) {
		return new GenomicIndexInterval( chromIdx, Math.max(1, start+offset), end+offset, getStrand() );
	}

	/**
	 * Create a new instance of GenomicInterval that is shifted from this one by "offset" (i.e. negative number to move upstream).
	 * @param offset
	 * @return shifted interval
	 */
	public GenomicIndexInterval createShiftedInterval(int offsetStart, int offsetEnd) {
		return new GenomicIndexInterval( chromIdx, Math.max(1, start+offsetStart), end+offsetEnd, getStrand() );
	}
}