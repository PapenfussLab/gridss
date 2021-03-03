package au.edu.wehi;

import java.io.File;

public interface VirusbreakenddbTests {
    static File findVirusbreakendFile(String file) {
        for (String rootPath : new String[] {
                "../virusbreakend/",
                "../virusbreakenddb_genbank_assemblies/",
                "../virusbreakenddb/",
        }) {
            for (String relativePath : new String[] {
                    "",
                    "taxonomy/",
                    "library/viral/",
                    "library/added/",
            }) {
                File f = new File(rootPath + relativePath + file);
                if (f.exists()) {
                    return f;
                }
            }
        }
        throw new RuntimeException("Cannot find reference genome to use for testing.");
    }
    static File getFullNodesDmp() { return findVirusbreakendFile("nodes.dmp"); }
    static File getSeqid2taxidMap() { return findVirusbreakendFile("seqid2taxid.map"); }
}
