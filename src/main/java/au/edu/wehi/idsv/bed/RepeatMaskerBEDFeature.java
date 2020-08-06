package au.edu.wehi.idsv.bed;

import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.bed.FullBEDFeature;

import java.awt.*;
import java.util.List;
import java.util.Map;

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

    public RepeatMaskerBEDFeature(Map<String, String> stringCache, BEDFeature feature, String[] tokens) {
        this.contig = getCached(stringCache, feature.getContig());
        this.start = feature.getStart();
        this.end = feature.getEnd();
        this.strand = feature.getStrand();
        this.swScore = feature.getScore();
        this.repeatType = getCached(stringCache, feature.getName());
        this.repeatClass = getCached(stringCache, tokens[10]);
    }

    private static String getCached(Map<String, String> stringCache, String s) {
        String cached = stringCache.get(s);
        if (cached == null) {
            cached = s;
            stringCache.put(s, s);
        }
        return cached;
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
