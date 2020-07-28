package au.edu.wehi.idsv.ncbi;

public class MinimalTaxonomyNode {
    public final int taxId;
    public final int parentTaxId;

    /**
     * Parse a line from nodes.dmp from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
     *
     * @param line line
     */
    public MinimalTaxonomyNode(String line) {
        int offset = 0;
        String[] fields = line.split("\t[|]\t");
        taxId = Integer.parseInt(fields[offset++]);
        parentTaxId = Integer.parseInt(fields[offset++]);
    }
}
