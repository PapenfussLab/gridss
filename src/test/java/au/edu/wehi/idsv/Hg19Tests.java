package au.edu.wehi.idsv;

import java.io.File;

/**
 * Tests requiring Hg19 
 * 
 * @author Daniel Cameron
 *
 */
public interface Hg19Tests {
	public static File findHg19Reference(String reference) {
		for (String path : new String[] {
				"examples/",
				"",
				"../",
				"C:/dev/",
				"~/projects/reference_genomes/human/",
			}) {
			File f = new File(path + reference);
			if (f.exists()) return f;
		}
		throw new RuntimeException("Cannot find hg19 reference genome to use for testing.");
	}
	public static File findHg19Reference() {
		return findHg19Reference("hg19.fa");
	}
}
