package gridss.cmdline;

import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;

import org.junit.Test;

import com.google.common.io.Files;

import au.edu.wehi.idsv.IntermediateFilesTest;

public class ReferenceCommandLineProgramTest extends IntermediateFilesTest {
	@Test
	public void should_create_dict() throws IOException {
		File fa = testFolder.newFile("test.fa");
		File dict = testFolder.newFile("test.fa.dict");
		dict.delete();
		Files.copy(SMALL_FA_FILE, fa);
		ReferenceCommandLineProgram.ensureSequenceDictionary(fa, null);
		assertTrue(dict.exists());
	}
}
