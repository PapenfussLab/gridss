package au.edu.wehi.idsv.bed;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.Iterator;

import org.apache.commons.io.LineIterator;

import com.google.common.collect.Iterators;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Log;

/**
 * Very basic parser that iterates over a BEDPE file
 * @author Daniel Cameron
 *
 */
public class BedpeIterator implements Iterator<BedpeRecord>, Closeable {
	private static final Log log = Log.getInstance(BedpeIterator.class);
	private BufferedReader br;
	private Iterator<BedpeRecord> it;
	public BedpeIterator(File file, SAMSequenceDictionary dict) throws FileNotFoundException {
		this(new FileReader(file), dict);
	}
	public BedpeIterator(Reader reader, SAMSequenceDictionary dict) {
		br = new BufferedReader(reader);
		Iterator<String> rawit = new LineIterator(br);
		final LineNumberContainer lineno = new LineNumberContainer();
		it = Iterators.transform(rawit, str -> {
			lineno.lineno++;
			// ignore headers and empty lines
			if (str == null || str.length() == 0 || str.charAt(0) == '#') {
				return null;
			}
			try {
				return new BedpeRecord(dict, str);
			} catch (IllegalArgumentException e) {
				log.info(String.format("Error parsing line %d. Ignoring", lineno.lineno));
				return null;
			} catch (Exception e) {
				log.error(e, String.format("Fatal error parsing record on line %d", lineno.lineno));
				throw e;
			}
		});
		it = Iterators.filter(it, br -> br != null);
	}
	private class LineNumberContainer {
		public int lineno = 0;
	}
	@Override
	public boolean hasNext() {
		return it.hasNext();
	}

	@Override
	public BedpeRecord next() {
		return it.next();
	}

	@Override
	public void close() throws IOException {
		br.close();
	}
}
