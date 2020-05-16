package au.edu.wehi.idsv.bed;

import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

public class RepeatMaskerBEDCodec extends BEDCodec {
    /**
     * Caching string to reduce memory usage of repeated strings
     */
    private Map<String, String> stringCache = new HashMap<>();
    @Override
    public BEDFeature decode(String[] tokens) {
        BEDFeature first6 = super.decode(Arrays.copyOf(tokens, 6));
        return new RepeatMaskerBEDFeature(stringCache, first6, tokens);
    }
}
