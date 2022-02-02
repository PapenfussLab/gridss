package assfolder;

import au.edu.wehi.idsv.IntermediateFilesTest;
import htsjdk.samtools.SAMRecord;
import org.junit.Test;

import java.io.File;
import java.util.List;

import static org.junit.Assert.assertTrue;

public class OlcExporterTest extends IntermediateFilesTest {
    @Test
    public void should_export_simple_overlap_graph() {
        File f = new File(testFolder.getRoot(),"test.OlcExporterTest.should_export_simple_overlap_graph.gexf");
        List<SAMRecord> reads = overlapping(1, 100, 20, 5);
        OlcExporter.exportOverlapGraph(reads, 5, 10, 0, f);
        assertTrue(f.exists());
        f.delete();
    }
    @Test
    public void should_export_gridss_example_graph() {
        File f = new File(testFolder.getRoot(),"test.OlcExporterTest.should_export_gridss_example_graph.gexf");
        List<SAMRecord> reads = getRecords(new File("src/test/resources/chr12.1527326.DEL1024.bam"));
        OlcExporter.exportOverlapGraph(reads, 20, 30, 2, f);
        assertTrue(f.exists());
    }
    /*
    @Test
    @Category(Hg19Tests.class)
    public void test_() {
        File f = new File("test.OlcExporterTest.colo829..gexf");
        SamReader reader = SamReaderFactory.makeDefault().open(new File("S:/colo829/COLO829T_dedup.realigned.bam"));
        List<SAMRecord> list = Lists.newArrayList();
        for (SAMRecord r : reader) {
            list.add(r);
        }
        List<SAMRecord> reads = getRecords(new File("example/chr12.1527326.DEL1024.bam"));
        OlcExporter.exportOverlapGraph(reads, 20, 30, 2, f);
        assertTrue(f.exists());
    }
    */
}