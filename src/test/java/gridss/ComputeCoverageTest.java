package gridss;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;

import org.junit.Assert;
import org.junit.Test;

import com.google.common.io.Files;

import au.edu.wehi.idsv.IntermediateFilesTest;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.readers.LineIterator;

public class ComputeCoverageTest extends IntermediateFilesTest {
	private static final double DELTA = 0.00001;
	private ArrayList<BEDFeature> getBed(File file) throws IOException {
		try (AbstractFeatureReader<BEDFeature, LineIterator> reader = AbstractFeatureReader.getFeatureReader(file.getAbsolutePath(), new BEDCodec(), false)) {
			Iterator<BEDFeature> it = reader.iterator();
			ArrayList<BEDFeature> list = new ArrayList<BEDFeature>();
			while (it.hasNext()) {
				list.add(it.next());
			}
			return list;
		}
	}
	private void expectBin(String chr, int start, int end, double value, BEDFeature bin) {
		Assert.assertEquals(chr, bin.getContig());
		Assert.assertEquals(start, bin.getStart());
		Assert.assertEquals(end, bin.getEnd());
		Assert.assertEquals(value, bin.getScore(), DELTA);
	}
	@Test
	public void should_return_average_bin_coverage() throws IOException {
		output = testFolder.newFile("out.bed");
		output.delete();
		createInput(RP(0, 5, 21, 5));
		String[] args = new String[] {
				"INPUT=" + input.toString(),
				"REFERENCE_SEQUENCE=" + reference.toString(),
				"OUTPUT=" + output.toString(),
				"TMP_DIR=" + super.testFolder.getRoot().toString(),
				"BIN_SIZE=10",
		};
		assertEquals(0, new ComputeCoverage().instanceMain(args));
		Assert.assertTrue(output.exists());
		// BEDFeature uses 1-based genomic coordinates 
		ArrayList<BEDFeature> list = getBed(output);
		expectBin("polyA", 1, 10, 0.5, list.get(0));
		expectBin("polyA", 11, 20, 0, list.get(1));
		expectBin("polyA", 21, 30, 0.5, list.get(2));
		expectBin("polyA", 31, 40, 0, list.get(3));
	}
	@Test
	public void should_adjust_For_gc_bias() throws IOException {
		output = new File(testFolder.getRoot(), "out.bed");
		File outputgc = new File(testFolder.getRoot(), "outgc.bed");
		File gcFile = new File(testFolder.getRoot(), "gcbias.txt");
		StringBuilder sb = new StringBuilder("0\t100\n"); // multiply 0 GC content by 100
		for (int i = 1; i <= 100; i++) {
			sb.append(i);
			sb.append('\t');
			sb.append(0.5);
			sb.append('\n');
		}
		Files.write(B(sb.toString()), gcFile);
		createInput(RP(0, 5, 21, 5));
		String[] args = new String[] {
				"INPUT=" + input.toString(),
				"REFERENCE_SEQUENCE=" + reference.toString(),
				"OUTPUT=" + output.toString(),
				"OUTPUT_GC=" + outputgc.toString(),
				"TMP_DIR=" + super.testFolder.getRoot().toString(),
				"BIN_SIZE=10",
				"GC_ADJUSTMENT=" + gcFile.toString(),
		};
		assertEquals(0, new ComputeCoverage().instanceMain(args));
		Assert.assertTrue(output.exists());
		// BEDFeature uses 1-based genomic coordinates 
		ArrayList<BEDFeature> list = getBed(output);
		expectBin("polyA", 1, 10, 0.5, list.get(0));
		expectBin("polyA", 11, 20, 0, list.get(1));
		expectBin("polyA", 21, 30, 0.5, list.get(2));
		expectBin("polyA", 31, 40, 0, list.get(3));
		
		Assert.assertTrue(outputgc.exists());
		// BEDFeature uses 1-based genomic coordinates 
		list = getBed(outputgc);
		expectBin("polyA", 1, 10, 50, list.get(0));
		expectBin("polyA", 11, 20, 0, list.get(1));
		expectBin("polyA", 21, 30, 50, list.get(2));
		expectBin("polyA", 31, 40, 0, list.get(3));
	}
}
