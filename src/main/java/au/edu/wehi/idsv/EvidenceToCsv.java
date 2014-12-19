package au.edu.wehi.idsv;

import htsjdk.samtools.util.CloseableIterator;

import java.io.BufferedOutputStream;
import java.io.Closeable;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.io.UnsupportedEncodingException;
import java.util.Iterator;

import com.google.common.collect.AbstractIterator;

/**
 * Writes the evidence stream to a file
 * @author Daniel Cameron
 *
 */
public class EvidenceToCsv extends AbstractIterator<DirectedEvidence> implements CloseableIterator<DirectedEvidence> {
	private final Iterator<DirectedEvidence> it;
	private final PrintStream stream;
	public EvidenceToCsv(File file, Iterator<DirectedEvidence> it) {
		this.it = it;
		try {
			this.stream = new PrintStream(new BufferedOutputStream(new FileOutputStream(file)));
		} catch (FileNotFoundException e) {
			throw new RuntimeException(e);
		}
		writeHeader();
	}
	private void writeHeader() {
	}
	@Override
	protected DirectedEvidence computeNext() {
		if (!it.hasNext()) return endOfData();
		DirectedEvidence e = it.next();
		double getBreakendQual();
		BreakendSummary getBreakendSummary();
		public byte[] getBreakendSequence();
		public byte[] getBreakendQuality();
		String getEvidenceID();
		EvidenceSource getEvidenceSource();
		int getLocalMapq();
		int getLocalBaseLength();
		int getLocalMaxBaseQual();
		int getLocalTotalBaseQual();
		
		stream.print("hello");
		return e;
	}
	@Override
	public void close() {
		stream.close();
	}
}
