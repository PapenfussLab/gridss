package au.edu.wehi.idsv;

import au.edu.wehi.idsv.configuration.GridssConfiguration;
import au.edu.wehi.idsv.picard.SynchronousReferenceLookupAdapter;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;

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
    static ProcessingContext createProcessingContext(File tempDir, File reference, GridssConfiguration config) {
        try {
            return new ProcessingContext(
                    new FileSystemContext(tempDir, tempDir, 500000),
                    reference,
                    new SynchronousReferenceLookupAdapter(new IndexedFastaSequenceFile(reference)), new ArrayList<>(),
                    config);
        } catch (FileNotFoundException e) {
            throw new RuntimeException(e);
        }
    }
}
