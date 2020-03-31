package au.edu.wehi.idsv;

import java.io.File;

/**
 * Tests requiring Hg19 
 * 
 * @author Daniel Cameron
 *
 */
public interface Hg38Tests extends ReferenceTests {
	static File findHg38Reference() {
		File f = ReferenceTests.findReference("Homo_sapiens_assembly38.fasta");
		if (f == null) {
			f = ReferenceTests.findReference("hg38.fa");
		}
		return f;
	}
}
