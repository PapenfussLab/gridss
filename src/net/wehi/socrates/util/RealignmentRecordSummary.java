/**
 * 
 */
package net.wehi.socrates.util;

import java.util.Arrays;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;

/**
 * @author hsu
 *
 * Created on Jan 24, 2013
 */
public class RealignmentRecordSummary {
	public int realignChr, realignPos, anchorChr, anchorPos, mapq;
	public boolean realignForward, anchorForward;
	public byte[] realignSeq, realignScSeq, anchorSeq;
	public boolean anchorIsIdeal;
//	public float alignedPercentIdentity;
	// realignSeq and anchorSeq have left most base adjacent to breakpoint, where realignScSeq has right most base adjacent to breakpoint
//	public String name;
	
	public RealignmentRecordSummary(SAMRecord realign) {
//		name = realign.getReadName().intern();
		// extract realignment info
		realignChr = realign.getReferenceIndex();
		byte[] seq = realign.getReadBases();
		mapq = realign.getMappingQuality();
		
		// determine if realigned sequence is soft-clipped at breakpoint
		java.util.List<CigarElement> cigar = realign.getCigar().getCigarElements();
		CigarElement breakpointCigar;
		
		if (realign.getReadNegativeStrandFlag()) {
			realignPos = realign.getAlignmentStart();
			realignForward = false;
			breakpointCigar = cigar.get(0);
		} else {
			realignPos = realign.getAlignmentEnd();
			Utilities.reverseComplement(seq);
			realignForward = true;
			breakpointCigar = cigar.get(cigar.size()-1);
		}
		
		int realignScBases = (breakpointCigar.getOperator() == CigarOperator.SOFT_CLIP) ? breakpointCigar.getLength() : 0;

		realignScSeq = Arrays.copyOfRange(seq, 0, realignScBases);
		realignSeq = Arrays.copyOfRange(seq, realignScBases, seq.length);
//		alignedPercentIdentity = calcAlignedPercentIdentity(realign);

		// reconstruct anchor info
		anchorForward = !(realign.getMateNegativeStrandFlag());
		anchorChr = realign.getMateReferenceIndex();
		anchorPos = realign.getMateAlignmentStart();
		anchorSeq = Utilities.stringToByte(realign.getStringAttribute("ZS"));
		anchorIsIdeal = !(realign.getMateUnmappedFlag());
	}
	
	public static String makeFastqHeader(SAMRecord aln, SAMRecordSummary summary, 
			int bpChr, int bpPos, boolean isIdeal, String alignSeq, char clipStrand) {
		StringBuilder name;
		if (aln.getReadPairedFlag())
			name = new StringBuilder(aln.getReadName() + (aln.getFirstOfPairFlag()?"/1":"/2"));
		else
			name = new StringBuilder(aln.getReadName());
		
        name.append( "&" + bpChr + "&" + bpPos );
//        name.append( "&" + (summary.isHeadClipped() ? "-" : "+") );
        name.append( "&" + clipStrand );
        name.append( "&" + (isIdeal?1:0) );
        name.append( "&" + alignSeq );
        return name.toString();
	}
	
	public String toString(SAMFileInfo header) {
		return (header.getSequenceName(realignChr) + ":" + realignPos + "\t" + (realignForward?'+':'-')) + "\t" + 
			   (header.getSequenceName(anchorChr) + ":" + anchorPos + "\t" + (anchorForward?'+':'-'));
	}

	public float calcAlignedPercentIdentity(SAMRecord aln) {
		Object nm = aln.getAttribute("NM");
		if (nm!=null) {
			return (float)(realignSeq.length-((Integer)nm).intValue())/(float)(realignSeq.length);
		}
		
    	Object md = aln.getAttribute("MD");
    	if (md==null) return 0;
    	
    	String MD = (String)md;
        int mismatch = 0;
        int match = 0;
        int tok = 0;
        StringBuilder b = new StringBuilder();
        for (int i=0; i<MD.length(); i++) {
            char c = MD.charAt(i);
            if (tok%2==0 && !Character.isDigit(c)) {
                tok++;
                match += Integer.parseInt( b.toString() );
                b = new StringBuilder();
            }
            if (tok%2==1 && Character.isDigit(c)) {
                tok++;
                mismatch += b.charAt(0)=='^' ? b.length()-1 : b.length();
                b = new StringBuilder();
            }
            b.append( c );
        }
        if (tok%2==0) match += Integer.parseInt( b.toString() );
        if (tok%2==1) mismatch += b.charAt(0)=='^' ? b.length()-1 : b.length();
        
        return (float)(match)/(match+mismatch);
    }
}
