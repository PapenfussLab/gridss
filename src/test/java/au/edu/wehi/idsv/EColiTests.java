package au.edu.wehi.idsv;

import java.io.File;

public interface EColiTests extends ReferenceTests {
    static File findEColiReference() {
        File f = ReferenceTests.findReference("BL21DE3_Genome.fasta");
        return f;
    }
}
