package au.edu.wehi.idsv.alignment;

import java.io.File;
import java.util.List;

import org.apache.commons.lang3.SystemUtils;

import gridss.SoftClipsToSplitReads;

/**
 * Tests requiring Hg19 
 * 
 * @author Daniel Cameron
 *
 */
public interface ExternalAlignerTests {
	public static final List<String> COMMAND_LINE = getDefaultAlignerCommandLine();
	public static final File REFERENCE = new File("chr12.fa");
	static List<String> getDefaultAlignerCommandLine() {
		List<String> cmd = new SoftClipsToSplitReads().ALIGNER_COMMAND_LINE;
		if (SystemUtils.IS_OS_WINDOWS) {
			// use WSL to invoke bwa on windows machines
			cmd.add(0, "wsl");
		}
		return cmd;
	}
}
