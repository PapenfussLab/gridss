package au.edu.wehi.idsv;

import au.edu.wehi.idsv.alignment.AlignerFactory;
import au.edu.wehi.idsv.alignment.SmithWatermanFastqAligner;

public class IterativeSplitReadRealignerTest extends SplitReadRealignerTest {
    private static final SmithWatermanFastqAligner aligner = new SmithWatermanFastqAligner(AlignerFactory.create(), 2);

    @Override
    protected SplitReadRealigner createAligner() {
        return new IterativeSplitReadRealigner(getContext(), aligner);
    }
}