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
		String seq = S(RANDOM);
		int readLength = 100;
		List<SAMRecord> reads = new ArrayList<>();
		for (int i = 0; i < seq.length(); i++) {
			String perfect_read = (seq + seq).substring(i, i + readLength);
			SAMRecord r = new SAMRecord(null);
			r.setReadName(Integer.toString(i));
			r.setReadBases(B(perfect_read));
			reads.add(r);
		}
		poc.exportOverlapGraph(reads, 10, output);
		Assert.assertTrue(output.exists());
	}
}
