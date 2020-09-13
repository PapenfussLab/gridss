package au.edu.wehi.idsv.util;

import com.google.common.io.Files;

import java.io.*;
import java.net.URI;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.List;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

public abstract class FileHelper {
	public static String[] GENOMIC_INDEX_FILES = new String[] {
		".bai",  // BAM
		".crai", // CRAM
		".idx",  // Tribble
		".tbi",  // tabix
		".csi",  // BCF
	};
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
		moveIndex(from, to, GENOMIC_INDEX_FILES);
	}
	public static void delete(File file, boolean deleteIndexes) throws IOException {
		file.delete();
		for (File f : getIndexFilesFor(file)) {
			f.delete();
		}
	}
	private static void moveIndex(File from, File to, String... indexSuffix) throws IOException {
		for (String index : indexSuffix) {
			moveIndex(from, to, index);
		}
	}
	private static void moveIndex(File from, File to, String indexSuffix) throws IOException {
		trymovesingle(
				new File(from.getPath() + indexSuffix),
				new File(to.getPath() + indexSuffix));
		
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
	public static void copy(File from, File to, boolean copyIndexes) throws IOException {
		if (!from.exists()) {
			throw new IllegalArgumentException("Cannot copy nonexist file" + from.getAbsolutePath());
		}
		if (to.exists()) {
			FileHelper.delete(to, copyIndexes);
		}
		Files.copy(from, to);
		copyIndex(from, to, GENOMIC_INDEX_FILES);
	}
	private static void copyIndex(File from, File to, String... indexSuffix) throws IOException {
		for (String index : indexSuffix) {
			copyIndex(from, to, index);
		}
	}
	private static void copyIndex(File from, File to, String indexSuffix) throws IOException {
		trycopysingle(
				new File(from.getPath() + indexSuffix),
				new File(to.getPath() + indexSuffix));
		
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
		List<File> index = new ArrayList<>();
		for (String indexSuffix : GENOMIC_INDEX_FILES) {
			File f = new File(file.getAbsolutePath() + indexSuffix);
			File f2 = new File(file.getParentFile(), Files.getNameWithoutExtension(file.getName()) + indexSuffix);
			if (f.exists()) index.add(f);
			if (f2.exists()) index.add(f2);
		}
		return index;
	}
	public static void zipDirectory(File zip, File directory) throws IOException {
		byte[] buffer = new byte[4194304];
		URI base = directory.toURI();
		Deque<File> queue = new ArrayDeque<>();
		queue.push(directory);
		try (OutputStream outStream = new FileOutputStream(zip)) {
			try (ZipOutputStream zipStream = new ZipOutputStream(outStream)) {
				while (!queue.isEmpty()) {
					directory = queue.pop();
					for (File child : directory.listFiles()) {
						String name = base.relativize(child.toURI()).getPath();
						if (child.isDirectory()) {
							queue.push(child);
							// Don't write directory entries
							//name = name.endsWith("/") ? name : name + "/";
							// zipStream.putNextEntry(new ZipEntry(name));
						} else {
							ZipEntry ze = new ZipEntry(name);
							ze.setSize(child.length());
							zipStream.putNextEntry(ze);
							try (FileInputStream fis = new FileInputStream(child)) {
								int len;
								while ((len = fis.read(buffer, 0, buffer.length)) > 0) {
									zipStream.write(buffer, 0, len);
								}
							}
							zipStream.closeEntry();
						}
					}
				}
			}
		}
	}
}
