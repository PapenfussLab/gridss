package au.edu.wehi.idsv;

import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import jaligner.Alignment;
import jaligner.Sequence;
import jaligner.SmithWatermanGotoh;
import jaligner.matrix.Matrix;

public class AlignmentHelper {
	private AlignmentHelper() { }
	static {
		// quieten down the JAligner logging spam
		java.util.logging.Logger.getLogger(SmithWatermanGotoh.class.getName());
		java.util.logging.LogManager.getLogManager().getLogger(SmithWatermanGotoh.class.getName()).setLevel(java.util.logging.Level.WARNING);
	}
	/**
	 * Converts an alignment of a read (seq2) against a reference sequence (seq1) to a read cigar
	 * @param alignment computed alignment
	 * @return CIGAR of seq2 relative to reference seq1
	 */
	public static Cigar alignmentToCigar(Alignment alignment) {
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
        	if (seq1[i] == Alignment.GAP) {
        		currentOp = CigarOperator.INSERTION;
        	} else if (seq2[i]  == Alignment.GAP) {
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
	public static Alignment align_local(Sequence ref, Sequence seq) {
		return align_bwamem_local(ref, seq);
	}
	/**
	 * Align using alignment algorithm similar to bowtie2
	 * @param ref reference sequence. Must be a small subsequence as full Smith Waterman Gotoh is performed 
	 * @param seq sequence 
	 * @return local sequence alignment
	 */
	public static Alignment align_bowtie2_local(Sequence ref, Sequence seq) {
		return SmithWatermanGotoh.align(ref, seq, bowtie2Matrix, bowtie2GapOpen, bowtie2GapExtend);
	}
	// bowtie2 scoring scheme for high-quality bases (http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#scores-higher-more-similar)
	private static final Matrix bowtie2Matrix = createMatrix(2, -6, -1);
	private static final float bowtie2GapOpen = 5;
	private static final float bowtie2GapExtend = 3;
	/**
	 * Align using alignment algorithm similar to bwa mem
	 * @param ref reference sequence. Must be a small subsequence as full Smith Waterman Gotoh is performed 
	 * @param seq sequence 
	 * @return local sequence alignment
	 */
	public static Alignment align_bwamem_local(Sequence ref, Sequence seq) {
		return SmithWatermanGotoh.align(ref, seq, bwamemMatrix, bwamemGapOpen, bwamemGapExtend);
	}
	private static final Matrix bwamemMatrix = createMatrix(1, -4, -4);
	private static final float bwamemGapOpen = 6;
	private static final float bwamemGapExtend = 1;
	
	private static Matrix createMatrix(int match, int mismatch, int ambiguous) {
		float[][] scores = new float[Matrix.SIZE][Matrix.SIZE];
        // Fill the matrix with the scores
        for (int i = 0; i < Matrix.SIZE; i++) {
            for (int j = 0; j < Matrix.SIZE; j++) {
                if (Character.toUpperCase(i) == Character.toUpperCase(j)) {
                    scores[i][j] = match;
                } else if (isUnambiguousBase(i) && isUnambiguousBase(j)) {
                    scores[i][j] = mismatch;
                } else {
                	scores[i][j] = ambiguous;
                }
            }
        }
        return new Matrix("bowtie2", scores);
	}
	private static boolean isUnambiguousBase(int base) {
		switch (base) {
		case 'A':
		case 'C':
		case 'G':
		case 'T':
			return false;
		default:
			return true;
		}
	}
}
