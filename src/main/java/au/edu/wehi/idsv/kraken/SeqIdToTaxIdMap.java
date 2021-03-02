package au.edu.wehi.idsv.kraken;

import htsjdk.samtools.util.RuntimeIOException;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public class SeqIdToTaxIdMap {
    public static Map<String, Integer> createLookup(File seqid2taxid_map) {
        try {
            List<String> lines = Files.readAllLines(seqid2taxid_map.toPath());
            return lines.stream()
                    .map(s -> s.split("\t"))
                    .collect(Collectors.toMap(s -> s[0], s -> Integer.parseInt(s[1])));
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }
}
