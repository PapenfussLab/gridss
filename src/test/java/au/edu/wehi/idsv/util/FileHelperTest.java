package au.edu.wehi.idsv.util;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;

import org.junit.Test;

import au.edu.wehi.idsv.IntermediateFilesTest;


public class FileHelperTest extends IntermediateFilesTest {
	@Test
	public void should_move() throws IOException {
		File from = testFolder.newFile();
		File to = testFolder.newFile();
		FileHelper.move(from, to, false);
		assertFalse(from.exists());
		assertTrue(to.exists());
	}
	@Test
	public void should_replace_destination() throws IOException {
		File from = testFolder.newFile();
		File to = new File(testFolder.getRoot(), "should_replace_destination.test");
		FileHelper.move(from, to, false);
		assertFalse(from.exists());
		assertTrue(to.exists());
	}
	@Test
	public void should_move_indicies() throws IOException {
		File from = new File(testFolder.getRoot(), "in.bam");
		File samtoolsbai = new File(testFolder.getRoot(), "in.bam.bai");
		File picardbai = new File(testFolder.getRoot(), "in.bai");
		File vcfidx = new File(testFolder.getRoot(), "in.bam.idx");
		File to = new File(testFolder.getRoot(), "out.bam");
		from.createNewFile();
		samtoolsbai.createNewFile();
		picardbai.createNewFile();
		vcfidx.createNewFile();
		FileHelper.move(from, to, true);
		assertFalse(from.exists());
		assertFalse(samtoolsbai.exists());
		assertFalse(picardbai.exists());
		assertFalse(vcfidx.exists());
		assertTrue(to.exists());
		assertTrue(new File(testFolder.getRoot(), "out.bai").exists());
		assertTrue(new File(testFolder.getRoot(), "out.bam.bai").exists());
		assertTrue(new File(testFolder.getRoot(), "out.bam.idx").exists());
	}
}