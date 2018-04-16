package htsjdk.samtools.fastq;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import com.google.common.collect.ImmutableList;

import au.edu.wehi.idsv.TestHelper;

public class NonFlushingBasicFastqWriterTest extends TestHelper {
	@Test
	public void should_match_BasicFastqWriter_output() throws IOException {
		TemporaryFolder testFolder = new TemporaryFolder();
		testFolder.create();
		List<FastqRecord> records = ImmutableList.of(
				new FastqRecord("read1", B("ACTG"), "qual1", B("1234")),
				new FastqRecord("r2", B("TGCAN"), "q2", B("12345"))
		);
		File bfq = testFolder.newFile("old.fq");
		File ufq = testFolder.newFile("new.fq");
		
		BasicFastqWriter bfw = new BasicFastqWriter(bfq);
		NonFlushingBasicFastqWriter fw = new NonFlushingBasicFastqWriter(ufq);
		for (FastqRecord fq : records) {
			bfw.write(fq);
			fw.write(fq);
		}
		bfw.close();
		fw.close();
		Assert.assertArrayEquals(Files.readAllBytes(bfq.toPath()), Files.readAllBytes(ufq.toPath()));
		testFolder.delete();
	}
}
