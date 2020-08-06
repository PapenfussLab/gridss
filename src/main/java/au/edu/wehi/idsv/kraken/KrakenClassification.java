package au.edu.wehi.idsv.kraken;

public class KrakenClassification {
    private final String line;
    public final boolean isClassified;
    public final String sequenceId;
    public final int taxonomyId;

    public String toKrakenOutput() { return line; }

    public KrakenClassification(String line) {
        this.line = line;
        String[] fields = line.split("\t");
        // "C"/"U": a one letter code indicating that the sequence was either classified or unclassified.
        this.isClassified = "C".equals(fields[0]);
        // The sequence ID, obtained from the FASTA/FASTQ header.
        this.sequenceId = fields[1];
        // The taxonomy ID Kraken 2 used to label the sequence; this is 0 if the sequence is unclassified.
        this.taxonomyId = Integer.parseInt(fields[2]);
        // The length of the sequence in bp. In the case of paired read data, this will be a string containing the lengths of the two sequences in bp, separated by a pipe character, e.g. "98|94".
        // Note that paired read data will contain a "|:|" token in this list to indicate the end of one read and the beginning of another.
        // fields[3] ignored
        // A space-delimited list indicating the LCA mapping of each k-mer in the sequence(s). For example, "562:13 561:4 A:31 0:1 562:3" would indicate that:
        // When Kraken 2 is run against a protein database (see [Translated Search]), the LCA hitlist will contain the results of querying all six frames of each sequence. Reading frame data is separated by a "-:-" token.
    }
    //public final IntList kmerRleLCA = new IntArrayList();
    //public final IntList kmerRleLength = new IntArrayList();
}
