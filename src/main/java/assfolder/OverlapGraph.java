package assfolder;

import java.util.Collection;

public class OverlapGraph {
    private final SeedLookup lookup;
    private final int maxMismatches;
    private final int minOverlap;

    public OverlapGraph(int seedSize, int minOverlap, int maxMismatches) {
        this.lookup = new SeedLookup(seedSize);
        if (minOverlap < seedSize) throw new IllegalArgumentException("minOverlap cannot be less than seedSize");
        this.minOverlap = minOverlap;
        this.maxMismatches = maxMismatches;
    }

    private void addSuccessor(Read left, Read right, int offset) {
        ReadOverlapSuccessor ros = new ReadOverlapSuccessor(left, right, offset);
        if (ros.mismatches > maxMismatches) return;
        if (ros.overlapLength < minOverlap) return;
        left.overlapSuccessors.add(ros);
    }

    /**
     * Adds the given read to the overlap graph
     * If the read already exists in the graph, the positional constraints
     * on the overlaps are relaxed
     * @param read Read to add
     * @return Updated or added read
     */
    public Read add(Read read) {
        Read existingRead = lookup.getRead(read);
        if (existingRead != null) {
            existingRead.relaxContstraints(read);
            read = existingRead;
            // TODO: findOverlaps restrictions
            if (read.constraints != null) throw new RuntimeException("NYI");
            //lookup.findAdditionalOverlaps(read, existingRead.constraints);
            return read;
        }
        Collection<ReadOffset> overlaps = lookup.findOverlaps(read);
        for (ReadOffset ro : overlaps) {
            if (ro.offset == 0) {
                addSuccessor(read, ro.read, 0);
                addSuccessor(ro.read, read, 0);
            } else if (ro.offset > 0) {
                addSuccessor(read, ro.read, ro.offset);
            } else {
                addSuccessor(ro.read, read, -ro.offset);
            }
        }
        lookup.add(read);
        return read;
    }
    public void remove(Read read) {
        // TODO: assert we are outside of the processing window
        lookup.remove(read);
    }
}
