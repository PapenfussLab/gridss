package au.edu.wehi.idsv.alignment;

import static org.junit.Assume.assumeTrue;

public class SswJniAlignerTest extends SmithWatermanAlignerTest {
    @Override
    protected Aligner create(int match, int mismatch, int ambiguous, int gapOpen, int gapExtend) {
        assumeTrue(AlignerFactory.isSswjniLoaded());
        return new SswJniAligner(match, mismatch, ambiguous, gapOpen, gapExtend);
    }
}