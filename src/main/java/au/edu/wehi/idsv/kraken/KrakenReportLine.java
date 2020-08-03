package au.edu.wehi.idsv.kraken;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.stream.Collectors;

public class KrakenReportLine {
    public double percentage;
    public long countAssignedToTree;
    public long countAssignedDirectly;
    public String rank;
    public int taxonomyId;
    public String scientificName;
    public KrakenReportLine(String line) {
        String[] fields = line.split("\t");
        int offset = 0;
        this.percentage = Double.parseDouble(fields[offset++].replace("%", ""));
        this.countAssignedToTree = Integer.parseInt(fields[offset++]);
        this.countAssignedDirectly = Integer.parseInt(fields[offset++]);
        this.rank = fields[offset++];
        this.taxonomyId = Integer.parseInt(fields[offset++]);
        this.scientificName = fields[offset++].trim();
    }
    public static List<KrakenReportLine> parseReport(File file) throws IOException {
        List<KrakenReportLine> report = Files.readAllLines(file.toPath()).stream()
                .map(s -> new KrakenReportLine(s))
                .collect(Collectors.toList());
        return report;
    }
}
