package au.edu.wehi.idsv.alignment;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.util.SequenceUtil;
import jaligner.Sequence;
import jaligner.SmithWatermanGotoh;
import jaligner.matrix.Matrix;

public class JAlignerAligner implements Aligner {
	private final Matrix matrix;
	private final int gapOpen;
	private final int gapExtend;
	public JAlignerAligner(int match, int mismatch, int ambiguous, int gapOpen, int gapExtend) {
		this.matrix = createMatrix(match, mismatch, ambiguous);
		this.gapOpen = gapOpen;
		this.gapExtend = gapExtend;
	}
	static {
		// quieten down the JAligner logging spam
		Logger jalignerLogger = Logger.getLogger(SmithWatermanGotoh.class.getName());
		jalignerLogger.setLevel(java.util.logging.Level.OFF);
	}
	@Override
	public Alignment align_smith_waterman(byte[] seq, byte[] ref) {
		Sequence sSeq = new Sequence("seq", new String(seq));
		Sequence sRef = new Sequence("ref", new String(ref));
		jaligner.Alignment jaln = SmithWatermanGotoh.align(sRef, sSeq, matrix, gapOpen, gapExtend);
		
		Alignment alignment = new Alignment(jaln.getStart1(), alignmentToCigar(jaln).toString());
		return alignment;
	}
	private static Matrix createMatrix(int match, int mismatch, int ambiguous) {
		float[][] scores = new float[Matrix.SIZE][Matrix.SIZE];
        // Fill the matrix with the scores
        for (int i = 0; i < Matrix.SIZE; i++) {
            for (int j = 0; j < Matrix.SIZE; j++) {
                if (Character.toUpperCase(i) == Character.toUpperCase(j)) {
                    scores[i][j] = match;
                } else if (SequenceUtil.isValidBase((byte) i) && SequenceUtil.isValidBase((byte) j)) {
                    scores[i][j] = mismatch;
                } else {
                	scores[i][j] = ambiguous;
                }
            }
        }
        return new Matrix("matrix", scores);
	}
	/**
	 * Converts an alignment of a read (seq2) against a reference sequence (seq1) to a read cigar
	 * @param alignment computed alignment
	 * @return CIGAR of seq2 relative to reference seq1
	 */
	private static Cigar alignmentToCigar(jaligner.Alignment alignment) {
		List<CigarElement> cigar = new ArrayList<CigarElement>();
        if (alignment.getStart2() > 0) {
        	cigar.add(new CigarElement(alignment.getStart2(), CigarOperator.SOFT_CLIP));
        }
        char[] seq1 = alignment.getSequence1();
        char[] seq2 = alignment.getSequence2();
        int length = 0;
        CigarOperator op = CigarOperator.MATCH_OR_MISMATCH;
        for (int i = 0; i < seq1.length; i++) {
        	CigarOperator currentOp;
        	if (seq1[i] == jaligner.Alignment.GAP) {
        		currentOp = CigarOperator.INSERTION;
        	} else if (seq2[i]  == jaligner.Alignment.GAP) {
        		currentOp = CigarOperator.DELETION;
        	} else {
        		currentOp = CigarOperator.MATCH_OR_MISMATCH;
        	}
        	if (currentOp != op) {
        		if (length > 0) {
        			cigar.add(new CigarElement(length, op));
        		}
        		op = currentOp;
        		length = 0;
        	}
        	length++;
        }
        cigar.add(new CigarElement(length, op));
        int basesConsumed = new Cigar(cigar).getReadLength();
        int seqLength = alignment.getOriginalSequence2().length();
        if (basesConsumed != seqLength) {
        	cigar.add(new CigarElement(seqLength - basesConsumed, CigarOperator.SOFT_CLIP));
        }
        return new Cigar(cigar);
	}
}
