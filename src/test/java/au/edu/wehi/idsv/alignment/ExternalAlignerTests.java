package au.edu.wehi.idsv.alignment;

import java.io.File;
import java.util.Arrays;
import java.util.List;

/**
 * Tests requiring Hg19 
 * 
 * @author Daniel Cameron
 *
 */
public interface ExternalAlignerTests {
	public static final List<String> COMMAND_LINE = Arrays.asList(new String[] {
			"C:/dev/bin/snap-beta.18-windows/snap.exe",
			"single",
			"C:/dev/bin/snap-beta.18-windows",
			"-wbs",
			"1",
			"-map",
			"-fastq",
			"-",
			"-o",
			"-sam",
			"-",
		});
	public static final File REFERENCE = new File("C:/dev/chr12.fa");
}
