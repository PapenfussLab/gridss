package au.edu.wehi.idsv.picard;

import au.edu.wehi.idsv.TestHelper;
import htsjdk.samtools.SAMSequenceRecord;
import org.junit.Assert;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import java.io.File;
import java.io.IOException;
import java.util.stream.Collectors;

import static org.junit.Assert.*;


public class TwoBitBufferedReferenceSequenceFileTest extends TestHelper {
	@Test
	public void getSequenceShouldMatchUnderlying() throws IOException {
		TwoBitBufferedReferenceSequenceFile b = new TwoBitBufferedReferenceSequenceFile(SMALL_FA);
		for (String contig : SMALL_FA.getSequenceDictionary().getSequences().stream().map(ssr -> ssr.getSequenceName()).collect(Collectors.toList())) {
			assertEquals(S(b.getSequence(contig).getBases()).toUpperCase(), S(SMALL_FA.getSequence(contig).getBases()).toUpperCase());
		}
		b.close();
	}
	@Test
	public void getSubsequenceAtShouldMatchUnderlying() throws IOException {
		TwoBitBufferedReferenceSequenceFile b = new TwoBitBufferedReferenceSequenceFile(SMALL_FA);
		for (String contig : SMALL_FA.getSequenceDictionary().getSequences().stream().map(ssr -> ssr.getSequenceName()).collect(Collectors.toList())) {
			for (int i = 1; i < 100; i++) {
				for (int j = i; j < 100; j++) {
					assertEquals(S(b.getSubsequenceAt(contig, i, j).getBases()).toUpperCase(), S(SMALL_FA.getSubsequenceAt(contig, i, j).getBases()).toUpperCase());
					assertEquals(b.getSubsequenceAt(contig, i, j).getName(), SMALL_FA.getSubsequenceAt(contig, i, j).getName());
					assertEquals(b.getSubsequenceAt(contig, i, j).getContigIndex(), SMALL_FA.getSubsequenceAt(contig, i, j).getContigIndex());
				}
			}
		}
		b.close();
	}
	@Test
	public void should_convert_ambiguous_bases_to_Ns() throws IOException {
		TwoBitBufferedReferenceSequenceFile b = new TwoBitBufferedReferenceSequenceFile(new InMemoryReferenceSequenceFile(new String[] { "test" }, new byte[][] { B("RYSWKMBDHVN.-") }));
		assertEquals("NNNNNNNNNNNNN", S(b.getSequence("test").getBases()));
		b.close();
	}
	@Test
	public void should_convert_to_uppercase() throws IOException {
		TwoBitBufferedReferenceSequenceFile b = new TwoBitBufferedReferenceSequenceFile(new InMemoryReferenceSequenceFile(new String[] { "test" }, new byte[][] { B("acgtn") }));
		assertEquals("ACGTN", S(b.getSequence("test").getBases()));
		b.close();
	}
	@Test
	public void should_lazily_cache() throws IOException {
		TemporaryFolder testFolder = new TemporaryFolder();
		testFolder.create();
		File file = new File(testFolder.getRoot(), "TwoBitBufferedReferenceSequenceFileTest.gridsscache");
		assertFalse(file.exists());
		TwoBitBufferedReferenceSequenceFile a = new TwoBitBufferedReferenceSequenceFile(SMALL_FA, file);
		assertFalse(file.exists());
		a.getBase(0, 1);
		assertTrue(file.exists());
	}
	@Test
	public void should_round_trip_through_cache() throws IOException {
		TemporaryFolder testFolder = new TemporaryFolder();
		testFolder.create();
		File file = new File(testFolder.getRoot(), "TwoBitBufferedReferenceSequenceFileTest.gridsscache");
		assertFalse(file.exists());
		TwoBitBufferedReferenceSequenceFile a = new TwoBitBufferedReferenceSequenceFile(SMALL_FA, file);
		a.getBase(0, 1);
		assertTrue(file.exists());
		TwoBitBufferedReferenceSequenceFile b = new TwoBitBufferedReferenceSequenceFile(SMALL_FA, file);
		for (int i = 0; i < SMALL_FA.getSequenceDictionary().size(); i++) {
			SAMSequenceRecord s = SMALL_FA.getSequenceDictionary().getSequence(i);
			assertEquals(
					S(SMALL_FA.getSequence(s.getSequenceName()).getBases()).toUpperCase(),
					S(b.getSequence(s.getSequenceName()).getBases()).toUpperCase());
			assertEquals(
					S(a.getSequence(s.getSequenceName()).getBases()),
					S(b.getSequence(s.getSequenceName()).getBases()));
		}
		file.delete();
		testFolder.delete();
	}
	@Test
	public void anyAmbiguous_should_return_ambiguous_overlap() throws IOException {
		TwoBitBufferedReferenceSequenceFile b = new TwoBitBufferedReferenceSequenceFile(new InMemoryReferenceSequenceFile(new String[] { "test" }, new byte[][] { B("ACGTNAAAAN") }));
		TwoBitBufferedReferenceSequenceFile.PackedReferenceSequence prs = b.getPackedSequence("test");
		Assert.assertFalse(prs.anyAmbiguous(1, 4));
		Assert.assertTrue(prs.anyAmbiguous(1, 4));
		Assert.assertTrue(prs.anyAmbiguous(4, 4));
		Assert.assertTrue(prs.anyAmbiguous(4, 5));
		Assert.assertTrue(prs.anyAmbiguous(3, 5));
		Assert.assertFalse(prs.anyAmbiguous(5, 5));
		Assert.assertTrue(prs.anyAmbiguous(5, 10));
	}
}
