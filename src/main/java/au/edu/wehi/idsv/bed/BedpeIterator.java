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

/**
 * Very basic parser that iterates over a BEDPE file
 * @author Daniel Cameron
 *
 */
public class BedpeIterator implements Iterator<BedpeRecord>, Closeable {
	private BufferedReader br;
	private Iterator<BedpeRecord> it;
	public BedpeIterator(File file, SAMSequenceDictionary dict) throws FileNotFoundException {
		this(new FileReader(file), dict);
	}
	public BedpeIterator(Reader reader, SAMSequenceDictionary dict) {
		br = new BufferedReader(reader);
		Iterator<String> rawit = new LineIterator(br);
		// strip headers and empty lines
		Iterator<String> filteredit = Iterators.filter(rawit, line -> line.length() > 0 && line.charAt(0) != '#');
		it = Iterators.transform(filteredit, str -> new BedpeRecord(dict, str));
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
