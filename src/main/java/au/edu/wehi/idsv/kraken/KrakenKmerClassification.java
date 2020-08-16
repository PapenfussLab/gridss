package au.edu.wehi.idsv.kraken;

import joptsimple.internal.Strings;

public class KrakenKmerClassification {
    public static final int AMBIGUOUS = -1;
    public final int taxonomyId;
    public final int kmerCount;
    public KrakenKmerClassification(String s) {
        String[] fields = s.split("[:]");
        taxonomyId = fields[0].equals("A") ? AMBIGUOUS : Integer.parseInt(fields[0]);
        kmerCount = (fields.length < 2 || Strings.isNullOrEmpty(fields[1])) ? 0 : Integer.parseInt(fields[1]);
    }
}
