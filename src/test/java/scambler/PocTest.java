package scambler;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;

import au.edu.wehi.idsv.IntermediateFilesTest;
import htsjdk.samtools.SAMRecord;

public class PocTest extends IntermediateFilesTest {
	@Test
	public void circular_example() {
		output = new File(testFolder.getRoot(), "uncompressed_string_graph.gexf");
		Poc poc = new Poc();
		String seq = S(RANDOM).substring(0, 100);
		String seq2 = S(RANDOM).substring(100, 200);
		int readLength = 50;
		List<SAMRecord> reads = new ArrayList<>();
		for (int i = 0; i < seq.length(); i += 1) {
			String perfect_read = (seq + seq2).substring(i, i + readLength);
			SAMRecord r = new SAMRecord(null);
			r.setReadName(Integer.toString(i));
			r.setReadBases(B(perfect_read));
			reads.add(r);
		}
		poc.exportOverlapGraph(reads, 16, output);
		Assert.assertTrue(output.exists());
	}
}
