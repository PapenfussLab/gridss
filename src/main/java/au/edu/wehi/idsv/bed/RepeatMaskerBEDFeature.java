package au.edu.wehi.idsv.bed;

import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.bed.FullBEDFeature;

import java.awt.*;
import java.util.Collections;
import java.util.List;

/**
 * See https://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/rmsk2bed.html#column-mapping for format
 */
public class RepeatMaskerBEDFeature implements BEDFeature {
    private final String contig;
    private final int start;
    private final int end;
    private final String repeatType;
    private final double swScore;
    private final Strand strand;
    private final String repeatClass;

    public RepeatMaskerBEDFeature(BEDFeature feature, String[] tokens) {
        this.start = feature.getStart();
        this.end = feature.getEnd();
        this.strand = feature.getStrand();
        this.swScore = feature.getScore();
        this.repeatType = feature.getName();
        this.contig = feature.getContig();
        this.repeatClass = tokens[10];
    }

    public String getRepeatType() {
        return repeatType;
    }

    public int getSmithWatermanScore() {
        return (int)swScore;
    }

    public String getRepeatClass() {
        return repeatClass;
    }

    @Override
    public Strand getStrand() {
        return strand;
    }

    @Override
    public String getType() {
        return getRepeatType();
    }

    @Override
    public Color getColor() { return null; }

    @Override
    public String getDescription() { return null; }

    @Override
    public List<FullBEDFeature.Exon> getExons()  { return null; }

    @Override
    public String getName() {
        return getRepeatType();
    }

    @Override
    public float getScore() {
        return getSmithWatermanScore();
    }

    @Override
    public String getLink() { return null; }

    @Override
    public String getContig() {
        return contig;
    }

    @Override
    public int getStart() {
        return start;
    }

    @Override
    public int getEnd() {
        return end;
    }
}
