package assfolder;

import au.edu.wehi.idsv.debruijn.PackedSequence;
import htsjdk.samtools.SAMRecord;

import java.util.Collection;
import java.util.HashSet;

public class Read extends PackedSequence {
    private final String readName; // do we need to hash this to reduce memory consumption?
    private final boolean onNegativeStrand;
    private final boolean secondInPair;
    public final PositionalConstraints constraints;
    public final Collection<ReadOverlapSuccessor> overlapSuccessors = new HashSet<>();
    public Read(SAMRecord record) {
        super(record.getReadBases(), false, false);
        this.readName = record.getReadName();
        this.onNegativeStrand = record.getReadNegativeStrandFlag();
        this.secondInPair = record.getReadPairedFlag() && record.getSecondOfPairFlag();
        this.constraints = null; // TODO
    }

    /**
     * Unique identified for this read.
     */
    public String uid() {
        // Potential optimisation: hash read name and encode in an int or long with 2 bits for strand/pairing
        return readName + (onNegativeStrand ? "-" : "") + (secondInPair ? "/2" : "/1");
    }

    public void relaxContstraints(Read read) {
        // TODO
    }

    @Override
    public boolean equals(Object obj) {
        if (!(obj instanceof Read)) return false;
        return uid().equals(((Read) obj).uid());
    }

    @Override
    public int hashCode() {
        return uid().hashCode();
    }
}
