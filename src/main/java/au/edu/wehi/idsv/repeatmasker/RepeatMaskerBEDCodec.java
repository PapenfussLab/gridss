package au.edu.wehi.idsv.repeatmasker;

import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * Parser for bedops RepeatMasker BED files
 * See https://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/rmsk2bed.html#column-mapping for format
 */
public class RepeatMaskerBEDCodec extends BEDCodec {
    private Map<String, String> stringCache = new HashMap<>();
    @Override
    public BEDFeature decode(String[] tokens) {
        BEDFeature first6 = super.decode(Arrays.copyOf(tokens, 6));
        RepeatMaskerFeature f = new RepeatMaskerFeature();
        f.setContig(first6.getContig());
        f.setStart(first6.getStart());
        f.setEnd(first6.getEnd());
        f.setStrand(first6.getStrand());
        f.setSwScore(first6.getScore());
        f.setRepeatType(first6.getName());
        if (tokens.length > 10) {
            f.setRepeatClass(tokens[10]);
        }
        f.useCachedString(stringCache);
        return f;
    }
}
