/**
 * 
 */
package au.edu.wehi.socrates.util;

import java.util.HashMap;
import java.util.LinkedList;

/**
 * @author hsu
 *
 * Created on Feb 7, 2013
 */
public class ConsensusMaker {
	private int matrixStart=-1, consensusPos=-1, consensusSupport=-1;
	private LinkedList<int[]> votingMatrix;
	private HashMap<Integer,Integer> breakpointCount;
	
	int increment;
	
	public ConsensusMaker(int breakPos, byte[] seq, boolean breakForward) {
		breakpointCount = new HashMap<Integer, Integer>();
		votingMatrix = new LinkedList<int[]>();
		if (seq.length==0) return;

		matrixStart = breakPos;
		
		breakpointCount.put( breakPos, 1 );
		
		for (int i=0; i<seq.length; i++) {
			int[] votes = new int[5];
			votes[ Utilities.byteToIndex(seq[i]) ] = 1;
			votingMatrix.add( votes );
		}
		
		increment = breakForward ? -1 : 1;
	}
	
	public void update(int breakPos, byte[] seq) {
		if (seq.length==0) return;
		
		// merge breakpoint position counts
		int count = breakpointCount.containsKey(breakPos) ? breakpointCount.get(breakPos)+1 : 1;
		breakpointCount.put( breakPos, count );
		if (matrixStart==-1) matrixStart = breakPos;
		
		int pos = breakPos;
		int csize = votingMatrix.size();
		LinkedList<int[]> prefixVotingMatrix = new LinkedList<int[]>();
		
		for (int i=0; i<seq.length; i++) {
			int baseIdx = Utilities.byteToIndex( seq[i] );
			int listIdx = (increment==-1) ? (matrixStart - pos) : (pos - matrixStart);
			if (listIdx < 0) {
				int[] vote = new int[5];
				vote[baseIdx] = 1;
				prefixVotingMatrix.add(vote);
			} else if (listIdx >= csize) {
				int[] vote = new int[5];
				vote[baseIdx] = 1;
				votingMatrix.addLast( vote );
			} else {
				votingMatrix.get( listIdx )[baseIdx] += 1;
			}
			pos += increment;
		}
		votingMatrix.addAll(0, prefixVotingMatrix);
		
		matrixStart = ((increment==-1) ? Math.max(breakPos, matrixStart) : Math.min(breakPos, matrixStart));
	}
	
	public int getConsensusBreakpointPosition() {
		if (consensusPos!=-1) return consensusPos;
		// consensus position
		int max_count = -1;
		int max_count_pos = -1;
		for (Integer p : breakpointCount.keySet()) {
			int count = breakpointCount.get(p);
			if (count > max_count) {
				max_count = count;
				max_count_pos = p;
			}
		}
		consensusPos = max_count_pos;
		consensusSupport = max_count;
		breakpointCount = null;
		
		return max_count_pos;
	}
	
	public int getConsensusSupport() {
		if (consensusSupport==-1) getConsensusBreakpointPosition();
		return consensusSupport;
	}
	
	public byte[] getConsensusSequence() {
		if (consensusPos!=-1) getConsensusBreakpointPosition();
		
		StringBuilder consensus = new StringBuilder();
		int pos = matrixStart;
		for (int[] vote : votingMatrix) {
			if (((increment==-1) && pos > consensusPos) || ((increment==1) && pos < consensusPos)) {
				pos += increment;
				continue; 
			}
			int max_vote = -1;
			int max_vote_idx = -1;
			int tot_vote = 0;
			for (int i=0; i<vote.length; i++) {
				tot_vote += vote[i];
				if (vote[i] > max_vote) {
					max_vote = vote[i];
					max_vote_idx = i;
				}
			}
			if ( (float)max_vote/(float)tot_vote <= 0.5 ) max_vote_idx = 0;
			char b = Utilities.indexToNucleotide(max_vote_idx);
			consensus.append(b);
			
			pos += increment;
		}
		votingMatrix = null;
		
		return Utilities.stringToByte(consensus.toString());
	}
	
	public void printConsensusMatrix() {
		if (votingMatrix==null) return;
		for (int i=0; i<5; i++) {
			System.out.print(Utilities.indexToNucleotide(i)+" ");
			for (int[] v : votingMatrix) {
				System.out.print(v[i] + " ");
			}
			System.out.println();
		}
	}
	
	
}
