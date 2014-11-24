package au.edu.wehi.idsv.visualisation;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.List;

import org.junit.Test;

import au.edu.wehi.idsv.IntermediateFilesTest;


public class SubgraphAssemblyAlgorithmTrackerBEDWriterTest extends IntermediateFilesTest {
	@Test
	public void should_sort_by_genomic_start_position() throws IOException {
		File f = new File(testFolder.getRoot(), "test.bed");
		SubgraphAssemblyAlgorithmTrackerBEDWriter writer = new SubgraphAssemblyAlgorithmTrackerBEDWriter(100, f);
		writer.write(new SubgraphAlgorithmMetrics(getContext(), 0, FWD) {{
			finalAnchors(12310, 12320);
		}});
		writer.write(new SubgraphAlgorithmMetrics(getContext(), 0, FWD) {{
			finalAnchors(12305, 12311);
		}});
		writer.write(new SubgraphAlgorithmMetrics(getContext(), 1, FWD) {{
			finalAnchors(12301, 12302);
		}});
		writer.close();
		List<String> list = Files.readAllLines(f.toPath(), StandardCharsets.US_ASCII);
		assertTrue(list.get(1).contains("12305"));
		assertTrue(list.get(2).contains("12310"));
		assertTrue(list.get(3).contains("12301"));
	}
	@Test
	public void should_write_track_header() throws IOException {
		File f = new File(testFolder.getRoot(), "test.bed");
		SubgraphAssemblyAlgorithmTrackerBEDWriter writer = new SubgraphAssemblyAlgorithmTrackerBEDWriter(100, f);
		writer.close();
		List<String> list = Files.readAllLines(f.toPath(), StandardCharsets.US_ASCII);
		assertEquals(1, list.size());
		assertTrue(list.get(0).startsWith("track"));
	}
}
