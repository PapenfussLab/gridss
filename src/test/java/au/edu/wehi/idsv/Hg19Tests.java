package au.edu.wehi.idsv;

import java.io.File;

/**
 * Tests requiring Hg19 
 * 
 * @author Daniel Cameron
 *
 */
public interface Hg19Tests extends ReferenceTests {
	static File findHg19Reference() {
		return ReferenceTests.findReference("hg19.fa");
	}
	static File findBroadHg19Reference() {
		return ReferenceTests.findReference("Homo_sapiens.GRCh37.GATK.illumina.fasta");
	}
}
