package au.edu.wehi.idsv.alignment;

public class JAlignerAlignerTest extends SmithWatermanAlignerTest {
    @Override
    protected Aligner create(int match, int mismatch, int ambiguous, int gapOpen, int gapExtend) {
        return new JAlignerAligner(match, mismatch, ambiguous, gapOpen, gapExtend);
    }
}