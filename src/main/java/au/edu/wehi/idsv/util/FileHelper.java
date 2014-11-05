package au.edu.wehi.idsv.util;

import java.io.File;
import java.io.IOException;

import com.google.common.io.Files;

public abstract class FileHelper {
	/**
	 * Moves the given file and any associated indexes 
	 * @param from
	 * @param to
	 * @throws IOException 
	 */
	public static void move(File from, File to, boolean moveIndexes) throws IOException {
		if (!from.exists()) {
			throw new IllegalArgumentException("Cannot move nonexist file" + from.getAbsolutePath());
		}
		if (to.exists()) {
			to.delete();
		}
		
		Files.move(from, to);
		moveIndex(from, to, ".bai");
		moveIndex(from, to, ".idx");
	}
	private static void moveIndex(File from, File to, String indexSuffix) throws IOException {
		trymovesingle(
				new File(from.getAbsolutePath() + indexSuffix),
				new File(to.getAbsolutePath() + indexSuffix));
		
		trymovesingle(new File(from.getParentFile(), Files.getNameWithoutExtension(from.getName()) + indexSuffix),
				new File(to.getParentFile(), Files.getNameWithoutExtension(to.getName()) + indexSuffix));
	}
	private static void trymovesingle(File from, File to) throws IOException {
		if (to.exists()) {
			to.delete();
		}
		if (from.exists()) {
			Files.move(from, to);
		}
	}
}
