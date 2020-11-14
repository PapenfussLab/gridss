package au.edu.wehi.idsv.ncbi;

import joptsimple.internal.Strings;

public enum TaxonomyLevel {
    //Unclassified("U", -1),
    Root("R", 0),
    Domain("D", 1),
    Kingdom("K", 2),
    Phylum("P", 3),
    Class("C", 4),
    Order("O", 5),
    Family("F", 6),
    Genus("G", 7),
    Species("S", 8);
    private final String krakenAbbreviation;
    private final int depth;
    TaxonomyLevel(String krakenAbbreviation, int depth) {
        this.krakenAbbreviation = krakenAbbreviation;
        this.depth = depth;
    }
    static public TaxonomyLevel parseKrakenReportRank(String rank) {
        if (Strings.isNullOrEmpty(rank)) return null;//Unclassified;
        rank = rank.substring(0, 1);
        for (TaxonomyLevel level : TaxonomyLevel.values()) {
            if (rank.equals(level.krakenAbbreviation)) {
                return level;
            }
        }
        throw new IllegalArgumentException(String.format("Unknown kraken rank %s", rank));
    }
    public int depth() { return depth; }
    public String krakenAbbreviation() { return krakenAbbreviation; }
}
