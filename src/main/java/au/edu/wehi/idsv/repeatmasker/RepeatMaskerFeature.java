package au.edu.wehi.idsv.repeatmasker;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.bed.FullBEDFeature;

import java.awt.*;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;

public class RepeatMaskerFeature implements BEDFeature {
    private String contig;
    private int start;
    private int end;
    private String repeatType;
    private float score;
    private Strand strand;
    private String repeatClass;
    private RepeatAlignmentSummaryInformation repeatAlignmentSummaryInformation;
    private String uid;
    private RepeatAlignmentInformation alignment;

    public RepeatMaskerFeature() {
    }

    /**
     * Reduces memory overhead of loading a RepeatMasker annotation file by
     * reusing common strings. Uncached strings are added to the cache lookup
     * except for the ID field as it is unique.
     *
     * Not required if performing streaming analysis
     * @param cache string cache
     */
    public void useCachedString(Map<String, String> cache) {
        this.contig = getCached(cache, contig);
        this.repeatType = getCached(cache, repeatType);
        this.repeatClass = getCached(cache, repeatClass);
        this.contig = getCached(cache, contig);
        if (cache.containsKey(uid)) {
            this.uid = cache.get(uid);
        }
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
    public float getScore() { return score; }

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

    public void setContig(String contig) {
        this.contig = contig;
    }

    public void setStart(int start) {
        this.start = start;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    public void setRepeatType(String repeatType) {
        this.repeatType = repeatType;
    }

    public void setScore(float score) {
        this.score = score;
    }

    public void setStrand(Strand strand) {
        this.strand = strand;
    }

    public void setRepeatClass(String repeatClass) {
        this.repeatClass = repeatClass;
    }

    public RepeatAlignmentSummaryInformation getRepeatAlignmentSummaryInformation() {
        return repeatAlignmentSummaryInformation;
    }

    public void setRepeatAlignmentSummaryInformation(RepeatAlignmentSummaryInformation rasi) {
        this.repeatAlignmentSummaryInformation = rasi;
    }

    public String getUniqueID() {
        return uid;
    }

    public void setUniqueID(String uid) {
        this.uid = uid;
    }

    public RepeatAlignmentInformation getRepeatAlignmentInformation(boolean infer) {
        if (alignment != null || !infer) return alignment;
        RepeatAlignmentInformation rai = new RepeatAlignmentInformation();
        rai.setRepeatStart(getRepeatAlignmentSummaryInformation().getMatchStart());
        List<CigarElement> ce = new ArrayList(2);
        if (getStart() > 1)  {
            ce.add(new CigarElement(getStart() - 1, CigarOperator.SOFT_CLIP));
        }
        ce.add(new CigarElement(getEnd() - getStart() + 1, CigarOperator.MATCH_OR_MISMATCH));
        rai.setCigar(new Cigar(ce));
        return rai;
    }

    public void setAlignment(RepeatAlignmentInformation alignment) {
        this.alignment = alignment;
    }

    public static final Comparator<RepeatMaskerFeature> ByAlignmentLength =
            Comparator.comparingInt((RepeatMaskerFeature f) -> Math.abs(f.getEnd() - f.getStart()));
    public static final Comparator<RepeatMaskerFeature> ByUniqueID = Comparator.comparing(RepeatMaskerFeature::getUniqueID);
}
