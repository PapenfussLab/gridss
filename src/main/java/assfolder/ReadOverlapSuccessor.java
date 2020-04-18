package assfolder;

import au.edu.wehi.idsv.debruijn.PackedSequence;

public class ReadOverlapSuccessor extends ReadOffset {
    public final int mismatches;
    public final int overlapLength;
    public final int trailingBases;
    public ReadOverlapSuccessor(Read left, Read right, int offset) {
        super(right, offset);
        assert(offset >= 0);
        this.overlapLength = PackedSequence.overlapLength(left, right, offset);
        this.mismatches = this.overlapLength - PackedSequence.overlapMatches(left, right, offset);
        this.trailingBases = right.length() - left.length() + offset;
    }
    public Read getSuccessor() {
        return read;
    }
    public int getLeadingBases() {
        return offset;
    }
    public int getTrailingBases() {
        return trailingBases;
    }
    public String getPrefixSequence(Read left) {
        return new String(left.getBytes(0, getLeadingBases()));
    }
    public String getSuffixSequence(Read left) {
        if (getTrailingBases() < 0) {
            return "-" + new String(left.getBytes(left.length() + getTrailingBases(), -getTrailingBases()));
        } else {
            return new String(read.getBytes(read.length() - getTrailingBases(), getTrailingBases()));
        }
    }
}
