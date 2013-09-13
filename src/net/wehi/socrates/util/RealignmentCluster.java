/**
 * 
 */
package net.wehi.socrates.util;

import java.util.Arrays;

/**
 * @author hsu
 *
 * Created on Jan 24, 2013
 */
public class RealignmentCluster implements Comparable<RealignmentCluster> {
	public int realignChr, realignPos, realignConsensusPos; 
	public int anchorChr, anchorPos, anchorConsensusPos;

	public boolean realignForward, anchorForward;
	public byte[] realignConsensusSCSeq = null, realignConsensusSeq = null, anchorConsensusSeq = null;
	public RealignmentCluster pairedCluster = null;
	
	public ConsensusMaker realignConsensusMaker, realignScConsensusMaker, anchorConsensusMaker;
	public int supportLong, supportShort, mapqTotal=0;
	public int supportLongBases=0, supportShortBases=0, supportShortMaxLen=0;
	public float mapqAvg = 0f;
	
	public RealignmentCluster(RealignmentRecordSummary realignment) {
		realignChr = realignment.realignChr;   realignForward = realignment.realignForward;
		anchorChr  = realignment.anchorChr;    anchorForward  = realignment.anchorForward;
		
		realignPos = realignment.realignPos; anchorPos = realignment.anchorPos;
		realignConsensusMaker = new ConsensusMaker(realignPos, realignment.realignSeq, realignForward);
		realignScConsensusMaker = new ConsensusMaker(realignForward?realignPos+1:realignPos-1, Utilities.getReversedComplementArray(realignment.realignScSeq), !realignForward);
		realignConsensusPos = -1;

		anchorConsensusMaker = new ConsensusMaker(anchorPos, realignment.anchorSeq, anchorForward);
		anchorConsensusPos = -1;
		
		supportLong = 1; supportShort = 0; mapqTotal += realignment.mapq;
		supportLongBases += realignment.realignSeq.length;
	}
	
	public RealignmentCluster(RealignmentRecordSummary realignment, int offset) {
		realignChr = realignment.realignChr;   realignForward = true;
		anchorChr  = realignment.anchorChr;    anchorForward  = true;
		realignPos = realignment.realignPos+offset; anchorPos = realignment.anchorPos+offset;
	}
	
	// make a reciprocal cluster
	public RealignmentCluster(RealignmentCluster cluster, int offset) {
		realignChr = cluster.anchorChr;   realignForward = !cluster.anchorForward; // flip orientation
		anchorChr  = cluster.realignChr;    anchorForward  = !cluster.realignForward;
		realignPos = cluster.anchorPos+offset; anchorPos = cluster.realignPos+offset;
	}
	
	// empty realignment cluster
	public RealignmentCluster(int realignChr, int realignPos, boolean realignFwd, int anchorChr, int anchorPos, boolean anchorFwd) {
		this.realignChr = realignChr; this.realignConsensusPos=realignPos; this.realignPos=realignPos; this.realignForward=realignFwd;
		this.anchorChr = anchorChr; this.anchorConsensusPos=anchorPos; this.anchorPos=anchorPos; this.anchorForward=anchorFwd;
		this.realignConsensusSeq = new byte[0]; this.realignConsensusSCSeq = new byte[0]; this.anchorConsensusSeq = new byte[0];
	}
	
	public void updateWith(RealignmentRecordSummary single) {
		// merge evidence counts
		supportLong++;
		mapqTotal += single.mapq;
		supportLongBases += single.realignSeq.length;

		realignConsensusMaker.update(single.realignPos, single.realignSeq);
		realignScConsensusMaker.update(realignForward?single.realignPos+1:single.realignPos-1, Utilities.getReversedComplementArray(single.realignScSeq));
		anchorConsensusMaker.update(single.anchorPos, single.anchorSeq);
	}
	
	public void callConsensusPosition() {
		// consensus position
		if (realignConsensusMaker==null || anchorConsensusMaker==null) return;
		realignConsensusPos = realignConsensusMaker.getConsensusBreakpointPosition();
		anchorConsensusPos = anchorConsensusMaker.getConsensusBreakpointPosition();
	}

	
	public void callConsensusSequence() {
		if (realignConsensusPos == -1 || anchorConsensusPos == -1) callConsensusPosition();
		if (realignConsensusSeq!=null || anchorConsensusSeq != null) return;
		
		realignConsensusSeq = realignConsensusMaker.getConsensusSequence();
		anchorConsensusSeq = anchorConsensusMaker.getConsensusSequence();
		
		int scpos = realignScConsensusMaker.getConsensusBreakpointPosition();
		if (scpos==-1) realignConsensusSCSeq = new byte[0];
		else {
			// number of soft clip needs to be at least half number of reads of consensus realign sequences
			int scSupport = realignScConsensusMaker.getConsensusSupport();
			int bpSupport = realignConsensusMaker.getConsensusSupport();
			if ((float)scSupport/(float)bpSupport < 0.5) {
				realignConsensusSCSeq = new byte[0];
				return;
			}
			
			int offset = realignForward ? (scpos - realignConsensusPos) : (realignConsensusPos - scpos);
			if (offset==1) {
				realignConsensusSCSeq = Utilities.getReversedComplementArray(realignScConsensusMaker.getConsensusSequence());
			} else {
				byte[] scs = Utilities.getReversedComplementArray(realignScConsensusMaker.getConsensusSequence());
				if (offset+1>=scs.length)
					realignConsensusSCSeq = new byte[0];
				else {
					int pos = (offset+1 >= 0) ? offset+1 : 0;
					realignConsensusSCSeq = Arrays.copyOfRange(scs, pos, scs.length);
				}
			}
		}
		
		realignConsensusMaker = null;
		anchorConsensusMaker = null;
	}
	
	public boolean nearTo(RealignmentRecordSummary realignment, int flank) {
		if (realignChr != realignment.realignChr || anchorChr != realignment.anchorChr) return false;
		if (realignForward != realignment.realignForward || anchorForward != realignment.anchorForward) return false;
		if ( Math.abs(realignPos - realignment.realignPos) > flank ) return false;
		if ( Math.abs(anchorPos - realignment.anchorPos) > flank ) return false;
		return true;
	}
	
	public boolean matchConsensus(RealignmentRecordSummary realignment, float similarity) {
		return true;
	}

	public boolean nearTo(RealignmentCluster realignment, int flank) {
		if (realignChr != realignment.realignChr || anchorChr != realignment.anchorChr) return false;
		if (realignForward != realignment.realignForward || anchorForward != realignment.anchorForward) return false;
		if ( Math.abs(realignPos - realignment.realignPos) > flank ) return false;
		if ( Math.abs(anchorPos - realignment.anchorPos) > flank ) return false;
		return true;
	}
	
	public boolean reciprocalNearTo(RealignmentCluster other, int flank) {
		if (realignChr != other.anchorChr || anchorChr != other.realignChr) return false;
		if (realignForward != other.anchorForward || anchorForward != other.realignForward) return false;
		if ( Math.abs(realignConsensusPos - other.anchorConsensusPos) > flank ) return false;
		if ( Math.abs(anchorConsensusPos - other.realignConsensusPos) > flank ) return false;
		return true;
	}
	
	@Override
	public int compareTo(RealignmentCluster other) {
		if (realignChr != other.realignChr) return realignChr - other.realignChr;
		if (realignPos != other.realignPos) return realignPos - other.realignPos;
		if (realignForward != other.realignForward) return (realignForward && !other.realignForward) ? -1 : 1;
		if (anchorChr != other.anchorChr) return anchorChr - other.anchorChr;
		if (anchorPos != other.anchorPos) return anchorPos - other.anchorPos;
		if (anchorForward != other.anchorForward) return (anchorForward && !other.anchorForward) ? -1 : 1;
		return 0;
	}
	
	@Override
	public boolean equals(Object o) {
		RealignmentCluster other = (RealignmentCluster)o;
		if (realignChr != other.realignChr) return false;
		if (realignConsensusPos != other.realignConsensusPos) return false;
		if (realignForward != other.realignForward) return false;
		if (anchorChr != other.anchorChr) return false;
		if (anchorConsensusPos != other.anchorConsensusPos) return false;
		if (anchorForward != other.anchorForward) return false;
		return true;
	}
	
	public String toString(SAMFileInfo info) {
		String rs = (realignConsensusSCSeq==null || realignConsensusSCSeq.length==0) ? Utilities.sequenceByteToString(realignConsensusSeq, false) :
			(Utilities.sequenceByteToString(realignConsensusSCSeq, false) + "*" + Utilities.sequenceByteToString(realignConsensusSeq, false));
		String as = Utilities.sequenceByteToString(anchorConsensusSeq, false);
		
		return info.getSequenceName(realignChr)+":"+realignConsensusPos+ "\t" + (realignForward?'+':'-') + "\t" + rs + "\t" +
			   info.getSequenceName(anchorChr) +":"+anchorConsensusPos + "\t" + (anchorForward? '+':'-') + "\t" + as + "\t" +
			   supportLong + "\t" + supportLongBases +"\t" + supportShort + "\t" + supportShortBases + "\t" + supportShortMaxLen + "\t" + mapqAvg;
	}

}
