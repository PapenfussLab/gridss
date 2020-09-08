package au.edu.wehi.idsv.repeatmasker;

import htsjdk.samtools.Cigar;

public class RepeatAlignmentInformation {
    private int repeatStart;
    private Cigar cigar;
    private int nestedBases;

    public int getRepeatStart() {
        return repeatStart;
    }

    public void setRepeatStart(int repeatStart) {
        this.repeatStart = repeatStart;
    }

    public Cigar getCigar() {
        return cigar;
    }

    public void setCigar(Cigar cigar) {
        this.cigar = cigar;
    }

    public void setNestedBases(int nestedBases) {
        this.nestedBases = nestedBases;
    }

    public int getNestedBases() {
        return nestedBases;
    }
}
