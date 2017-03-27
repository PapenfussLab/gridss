package au.edu.wehi.idsv;

import java.io.File;

/**
 * Tests requiring Hg19 
 * 
 * @author Daniel Cameron
 *
 */
public interface Hg38Tests {
	public static File findHg38Reference(String reference) {
		for (String path : new String[] {
				"../",
				"C:/dev/",
				"C:/dev/",
				"~/projects/reference_genomes/human/",
			}) {
			File f = new File(path + reference);
			if (f.exists()) return f;
		}
		throw new RuntimeException("Cannot find hg38 reference genome to use for testing.");
	}
	public static File findHg38Reference() {
		return findHg38Reference("Homo_sapiens_assembly38.fasta");
	}
}
