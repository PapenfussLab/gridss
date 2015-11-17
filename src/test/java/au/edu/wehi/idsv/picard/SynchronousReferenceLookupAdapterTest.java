package au.edu.wehi.idsv.picard;

import static org.junit.Assert.assertEquals;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

import java.io.File;
import java.util.Random;

import org.junit.Test;

public class SynchronousReferenceLookupAdapterTest {
	private Exception e;
	private File findBigReference() {
		for (String path : new String[] {
				"C:/dev/hg19_karyotypic.fa",
				"C:/dev/hg19.fa",
				"~/projects/reference_genomes/human/hg19.fa",
			}) {
			File f = new File(path);
			if (f.exists()) return f;
		}
		throw new RuntimeException("Cannot file large reference genome to use for testing.");
	}
	@Test
	public void multi_threaded_read() throws Exception {
		File ref = findBigReference();
		IndexedFastaSequenceFile indexed = new IndexedFastaSequenceFile(ref);
		SynchronousReferenceLookupAdapter a = new SynchronousReferenceLookupAdapter(indexed);
		e = null;
		Thread[] th = new Thread[16];
		for (int t = 0; t < th.length; t++) {
			th[t] = new Thread(() -> {
				try {
					Random rng = new Random();
					for (int i = 0; i < 4096; i++) {
						int pos = rng.nextInt(100000000);
						assertEquals(a.getBase(indexed.getSequenceDictionary().getSequence("chr4").getSequenceIndex(),pos),
								indexed.getSubsequenceAt("chr4", pos, pos).getBases()[0]);
					}
				} catch (Exception ex) {
					e = ex;
				}
			});
			th[t].start();
		}
		for (int t = 0; t < th.length; t++) {
			th[t].join();
		}
		if (e != null) throw e;
		a.close();
	}
}
