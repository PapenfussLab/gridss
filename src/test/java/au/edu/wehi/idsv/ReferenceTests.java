package au.edu.wehi.idsv;

import java.io.File;

public interface ReferenceTests {
    static File findReference(String reference) {
        for (String path : new String[] {
                "examples/",
                "",
                "../",
                "../../",
                "ref/",
                "../ref/",
                "../../ref/",
                "C:/dev/",
                "D:/dev/",
                "~/projects/reference_genomes/human/",
        }) {
            File f = new File(path + reference);
            if (f.exists()) {
                return f;
            }
        }
        throw new RuntimeException("Cannot find reference genome to use for testing.");
    }
}
