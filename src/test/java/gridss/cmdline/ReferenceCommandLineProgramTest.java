package gridss.cmdline;

import au.edu.wehi.idsv.IntermediateFilesTest;
import com.google.common.io.Files;
import org.junit.Test;

import java.io.File;
import java.io.IOException;

import static org.junit.Assert.assertTrue;

public class ReferenceCommandLineProgramTest extends IntermediateFilesTest {
	@Test
	public void should_create_dict() throws IOException {
		File fa = testFolder.newFile("test.fa");
		File dict = testFolder.newFile("test.fa.dict");
		dict.delete();
		Files.copy(SMALL_FA_FILE, fa);
		ReferenceCommandLineProgram.ensureSequenceDictionary(fa);
		assertTrue(dict.exists());
	}
}
