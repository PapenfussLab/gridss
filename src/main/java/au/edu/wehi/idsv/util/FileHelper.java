package au.edu.wehi.idsv.util;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

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
			FileHelper.delete(to, moveIndexes);
		}
		if (!from.renameTo(to)) {
			throw new IOException("Could not rename " + from + " to " + to);
		}
		moveIndex(from, to, ".bai");
		moveIndex(from, to, ".idx");
	}
	public static void delete(File file, boolean deleteIndexes) throws IOException {
		file.delete();
		for (File f : getIndexFilesFor(file)) {
			f.delete();
		}
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
	/**
	 * Moves the given file and any associated indexes 
	 * @param from
	 * @param to
	 * @throws IOException 
	 */
	public static void copy(File from, File to, boolean moveIndexes) throws IOException {
		if (!from.exists()) {
			throw new IllegalArgumentException("Cannot copy nonexist file" + from.getAbsolutePath());
		}
		if (to.exists()) {
			FileHelper.delete(to, moveIndexes);
		}
		Files.copy(from, to);
		copyIndex(from, to, ".bai");
		copyIndex(from, to, ".idx");
	}
	private static void copyIndex(File from, File to, String indexSuffix) throws IOException {
		trycopysingle(
				new File(from.getAbsolutePath() + indexSuffix),
				new File(to.getAbsolutePath() + indexSuffix));
		
		trycopysingle(new File(from.getParentFile(), Files.getNameWithoutExtension(from.getName()) + indexSuffix),
				new File(to.getParentFile(), Files.getNameWithoutExtension(to.getName()) + indexSuffix));
	}
	private static void trycopysingle(File from, File to) throws IOException {
		if (to.exists()) {
			to.delete();
		}
		if (from.exists()) {
			Files.copy(from, to);
		}
	}
	public static List<File> getIndexFilesFor(File file) {
		return Stream.concat(
				getPossibleIndexFilesFor(file, ".bai"),
				getPossibleIndexFilesFor(file, ".idx"))
			.filter(f -> f.exists())
			.collect(Collectors.toList());
	}
	private static Stream<File> getPossibleIndexFilesFor(File file, String indexSuffix) {
		return Stream.of(
				new File(file.getAbsolutePath() + indexSuffix),
				new File(file.getParentFile(), Files.getNameWithoutExtension(file.getName()) + indexSuffix));
	}
}
