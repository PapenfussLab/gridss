package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;

import java.io.File;

import org.junit.Test;


public class EvidenceSourceTest extends IntermediateFilesTest {
	public class TestEvidenceSource extends EvidenceSource {
		public TestEvidenceSource(boolean perChr, File input) {
			super(getCommandlineContext(perChr), input);
		}
		protected boolean processed = false;
		@Override
		public int getMaxConcordantFragmentSize() {
			return 300;
		}
	}
	@Test
	public void getFileIntermediateDirectoryBasedOn_should_be_input_file() {
		EvidenceSource es = new TestEvidenceSource(false, input);
		assertEquals(input, es.getFileIntermediateDirectoryBasedOn());
	}
	
}
