package au.edu.wehi.idsv.bed;

import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.bed.FullBEDFeature;

import java.awt.*;
import java.util.List;

public class RepeatMaskerBEDFeature implements BEDFeature {
    private final BEDFeature feature;
    private final String[] tokens;

    public String getRepeatType() {
        return getName();
    }

    public int getSmithWatermanScore() {
        return (int)getScore();
    }

    public String getRepeatClass() {
        return tokens[10];
    }

    public RepeatMaskerBEDFeature(BEDFeature feature, String[] tokens) {
        this.feature = feature;
        this.tokens = tokens;
    }

    @Override
    public Strand getStrand() {
        return feature.getStrand();
    }

    @Override
    public String getType() {
        return feature.getType();
    }

    @Override
    public Color getColor() {
        return feature.getColor();
    }

    @Override
    public String getDescription() {
        return feature.getDescription();
    }

    @Override
    public List<FullBEDFeature.Exon> getExons() {
        return feature.getExons();
    }

    @Override
    public String getName() {
        return feature.getName();
    }

    @Override
    public float getScore() {
        return feature.getScore();
    }

    @Override
    public String getLink() {
        return feature.getLink();
    }

    @Override
    public String getContig() {
        return feature.getContig();
    }

    @Override
    public int getStart() {
        return feature.getStart();
    }

    @Override
    public int getEnd() {
        return feature.getEnd();
    }
}
