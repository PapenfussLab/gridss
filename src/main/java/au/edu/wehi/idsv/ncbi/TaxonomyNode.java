package au.edu.wehi.idsv.ncbi;

public class TaxonomyNode extends MinimalTaxonomyNode {
    public final String rank;
    public final String emblCode;
    public final int divisionId;
    public final boolean inheritedDivision;
    public final int geneticCodeID;
    public final boolean inheritedGeneticCodeID;
    public final int mitochondrialGeneticCodeId;
    public final boolean inheritedMitochondrialGeneticCodeId;
    public final boolean hiddenInGenBank;
    public final boolean hiddenSubTreeRoot;
    public final String comments;

    /**
     * Parse a line from nodes.dmp from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
     *
     * @param line line
     */
    public TaxonomyNode(String line) {
        super(line);
        int offset = 2;
        String[] fields = line.split("\t[|]\t");
        rank = fields[offset++];
        emblCode = fields[offset++];
        divisionId = Integer.parseInt(fields[offset++]);
        inheritedDivision = Integer.parseInt(fields[offset++]) == 1;
        geneticCodeID = Integer.parseInt(fields[offset++]);
        inheritedGeneticCodeID = Integer.parseInt(fields[offset++]) == 1;
        mitochondrialGeneticCodeId = Integer.parseInt(fields[offset++]);
        inheritedMitochondrialGeneticCodeId = Integer.parseInt(fields[offset++]) == 1;
        hiddenInGenBank = Integer.parseInt(fields[offset++]) == 1;
        hiddenSubTreeRoot = Integer.parseInt(fields[offset++]) == 1;
        comments = stripTrailingTabColon(fields[offset++]);
    }
    private static String stripTrailingTabColon(String s) {
        return s.substring(0, s.length() - 2);
    }
}
