package au.edu.wehi.idsv.kraken;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

public class KrakenReportLine {
    public final String line;
    public final double percentage;
    public final long countAssignedToTree;
    public final long countAssignedDirectly;
    public final String rank;
    public final int taxonomyId;
    public final String scientificName;
    public KrakenReportLine(String line) {
        this.line = line;
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
    public static final Comparator<KrakenReportLine> ByCountAssignedDirectly = (o1, o2) -> Long.compare(o1.countAssignedDirectly, o2.countAssignedDirectly);
    public static final Comparator<KrakenReportLine> ByCountAssignedToTree = (o1, o2) -> Long.compare(o1.countAssignedToTree, o2.countAssignedToTree);
}
