package au.edu.wehi.idsv.picard;

import au.edu.wehi.idsv.Hg19Tests;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import org.junit.Test;
import org.junit.experimental.categories.Category;

import java.io.File;
import java.util.Random;

import static org.junit.Assert.assertEquals;

public class SynchronousReferenceLookupAdapterTest {
	private Exception e;
	@Test
	@Category(Hg19Tests.class)
	public void multi_threaded_read() throws Exception {
		File ref = Hg19Tests.findHg19Reference();
		IndexedFastaSequenceFile indexed = new IndexedFastaSequenceFile(ref);
		SynchronousReferenceLookupAdapter a = new SynchronousReferenceLookupAdapter(indexed);
		e = null;
		Thread[] th = new Thread[16];
		for (int t = 0; t < th.length; t++) {
			th[t] = new Thread(() -> {
				try {
					Random rng = new Random(0);
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
		a.close();
		if (e != null) throw e;
	}
}
