package au.edu.wehi.idsv.debruijn.positional;

import com.google.common.collect.Range;

/**
 * Assembly graph has become too dense in the given region
 */
public class AssemblyThresholdReachedException extends RuntimeException {
    private final Range<Integer> range;

    public AssemblyThresholdReachedException(Range<Integer> range) {
        this.range = range;
    }

    public Range<Integer> getRange() {
        return range;
    }
}
