package au.edu.wehi.idsv.bed;

import java.io.BufferedOutputStream;
import java.io.Closeable;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;

import htsjdk.samtools.SAMSequenceDictionary;

/**
 * Writes to a BED file
 * @author Daniel Cameron
 *
 */
public class BedWriter implements Closeable {
	private OutputStream os;
	private SAMSequenceDictionary dict;
	public BedWriter(SAMSequenceDictionary dictionary, File file) throws IOException {
		this(dictionary, new FileOutputStream(file));
	}
	public BedWriter(SAMSequenceDictionary dictionary, OutputStream stream) throws IOException {
		this.dict = dictionary;
		this.os = new BufferedOutputStream(stream);
	}
	public void writeHeader() throws IOException {
	}
	/**
	 * Writes the given information to a BED file 
	 * @param referenceIndex reference index
	 * @param start 1-based inclusive start position of interval 
	 * @param end 1-based inclusive end position of interval
	 * @param score score
	 * @throws IOException
	 */
	public void write(int referenceIndex, int start, int end, double score) throws IOException {
		StringBuilder sb = new StringBuilder();
		sb.append(dict.getSequence(referenceIndex).getSequenceName());
		sb.append('\t'); sb.append(Integer.toString(start - 1));
		sb.append('\t'); sb.append(Integer.toString(end - 1));
		sb.append('\t'); sb.append(Double.toString(score));
		sb.append('\n');
		os.write(sb.toString().getBytes(StandardCharsets.UTF_8));
	}
	@Override
	public void close() throws IOException {
		if (os != null) os.close();
		os = null;
	}
}
