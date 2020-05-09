package au.edu.wehi.idsv.bed;

import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;

import java.util.Arrays;

public class RepeatMaskerBEDCodec extends BEDCodec {
    @Override
    public BEDFeature decode(String[] tokens) {
        BEDFeature first6 = super.decode(Arrays.copyOf(tokens, 6));
        return new RepeatMaskerBEDFeature(first6, tokens);
    }
}
