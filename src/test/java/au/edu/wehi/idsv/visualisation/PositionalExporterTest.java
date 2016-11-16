package au.edu.wehi.idsv.visualisation;

import static org.junit.Assert.assertTrue;

import java.io.IOException;

import org.junit.Test;

import com.google.common.collect.ImmutableList;

import au.edu.wehi.idsv.IntermediateFilesTest;
import au.edu.wehi.idsv.debruijn.positional.KmerPathNode;


public class PositionalExporterTest extends IntermediateFilesTest {
	@Test
	public void should_export_single_node() throws IOException {
		PositionalExporter.exportDot(output, 4, ImmutableList.of(new KmerPathNode(0, 1, 2, false, 3)), null);
		assertTrue(output.exists());
	}
	@Test
	public void should_export_full_node() throws IOException {
		PositionalExporter.exportNodeDot(output, 4, ImmutableList.of(new KmerPathNode(0, 1, 2, false, 3)), null);
		assertTrue(output.exists());
	}
}
