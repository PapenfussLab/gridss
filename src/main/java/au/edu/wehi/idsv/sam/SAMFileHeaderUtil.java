package au.edu.wehi.idsv.sam;

import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMFileHeader;

public class SAMFileHeaderUtil {
    /**
     * Minimal header required so merging tools don't duplicate headers such as comments
     */
    public static SAMFileHeader minimal(SAMFileHeader header) {
        header = header.clone();
        // Strip everything that could be duplicated in the merger
        header.setComments(ImmutableList.of());
        header.setProgramRecords(ImmutableList.of());
        return header;
    }
}
