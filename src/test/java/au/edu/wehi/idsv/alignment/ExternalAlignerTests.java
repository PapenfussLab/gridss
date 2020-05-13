package au.edu.wehi.idsv.alignment;

import au.edu.wehi.idsv.Hg19Tests;
import au.edu.wehi.idsv.ReferenceTests;
import gridss.SoftClipsToSplitReads;
import org.apache.commons.lang3.SystemUtils;

import java.io.File;
import java.util.List;

/**
 * Tests requiring Hg19 
 * 
 * @author Daniel Cameron
 *
 */
public interface ExternalAlignerTests {
	public static final List<String> COMMAND_LINE = getDefaultAlignerCommandLine();
	public static final File REFERENCE = ReferenceTests.findReference("chr12.fa");
	static List<String> getDefaultAlignerCommandLine() {
		List<String> cmd = new SoftClipsToSplitReads().ALIGNER_COMMAND_LINE;
		if (SystemUtils.IS_OS_WINDOWS) {
			// use WSL to invoke bwa on windows machines
			cmd.add(0, "wsl");
		}
		return cmd;
	}
}
