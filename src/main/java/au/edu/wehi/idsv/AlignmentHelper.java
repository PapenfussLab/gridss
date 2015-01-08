package au.edu.wehi.idsv;

import jaligner.Alignment;
import jaligner.Sequence;
import jaligner.SmithWatermanGotoh;
import jaligner.matrix.Matrix;

public class AlignmentHelper {
	private AlignmentHelper() { }
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
                if (i == j) {
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
