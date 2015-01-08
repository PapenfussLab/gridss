package au.edu.wehi.idsv;

import jaligner.Alignment;
import jaligner.Sequence;
import jaligner.SmithWatermanGotoh;
import jaligner.matrix.Matrix;

public class AlignmentHelper {
	private AlignmentHelper() { }
	/**
	 * Align using alignment algorithm similar to bowtie2
	 * @param ref reference sequence. Must be a small subsequence as full Smith Waterman Gotoh is performed 
	 * @param seq sequence 
	 * @return local sequence alignment
	 */
	public static Alignment align_bowtie2_local(Sequence ref, Sequence seq) {
		return SmithWatermanGotoh.align(ref, seq, bowtie2Matrix, bowtie2GapOpen, bowtie2GapExtend);
	}
	// Default to bowtie2 scoring scheme for high-quality bases (http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#scores-higher-more-similar)
	private static final Matrix bowtie2Matrix = createMatrixBowtie2();
	private static final float bowtie2GapOpen = 5;
	private static final float bowtie2GapExtend = 3;
	private static Matrix createMatrixBowtie2() {
		float[][] scores = new float[Matrix.SIZE][Matrix.SIZE];
        // Fill the matrix with the scores
        for (int i = 0; i < Matrix.SIZE; i++) {
            for (int j = 0; j < Matrix.SIZE; j++) {
                if (i == j) {
                    scores[i][j] = 1; // 2
                } else if (isUnambiguousBase(i) && isUnambiguousBase(j)) {
                    scores[i][j] = -6;
                } else {
                	scores[i][j] = -1;
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
