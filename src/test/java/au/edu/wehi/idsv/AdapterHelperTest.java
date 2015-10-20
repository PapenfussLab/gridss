package au.edu.wehi.idsv;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import htsjdk.samtools.SAMRecord;

import org.junit.Test;

import au.edu.wehi.idsv.configuration.SoftClipConfiguration;


public class AdapterHelperTest extends TestHelper {
	@Test
	public void containsAdapter() {
		AdapterHelper ah = new SoftClipConfiguration().adapters;
		SAMRecord r = Read(0,  1,  101);
		r.setReadBases(B("CATGCCTGGCGTCCACTTTCTTGACTATTTCCTGAAGACCAGCGTTTCCCGGGTGGGTTCACAGCTGCGGAAGCTGCCTGTGTCAAGATCGGAAGAGCGTC"));
		assertTrue(ah.containsAdapter(r));
		
		r.setReadBases(B("ACCGCTCTTCCGATCTCATGCCTGGCGTCCACTTTCTTGACTATTTCCTGAAGACCAGCGTTTCCCGGGTGGTTTCACAGCTGCGGAAGCTGCCTGTGTCA"));
		assertTrue(ah.containsAdapter(r));
		
		r.setReadBases(B("TGAGCCACCATGCCTGGCGTCCACTTTCTTGACTATTTCCTGAAGACCAGCGTTTCCCGGGTGGTTTCACAGCTGCGGAAGCTGCCTGTGTCAGTGTCAAG"));
		assertFalse(ah.containsAdapter(r));
	}
}
