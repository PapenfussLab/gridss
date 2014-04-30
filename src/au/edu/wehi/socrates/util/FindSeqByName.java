/**
 * 
 */
package au.edu.wehi.socrates.util;

import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

/**
 * @author hsu
 *
 * Created on Feb 28, 2013
 */
public class FindSeqByName implements Callable<SAMRecord> {
	private String filename = null;
	private String readId = null;
	private String chrom = null;
	private int chromLen = 0;
	
	public static void findSeqByName(String bamFilename, String readName) {
		findSeqByName(bamFilename, readName, 3);
	}
	
	public static void findSeqByName(String bamFilename, String readName, int threads) {
		ExecutorService pool = null;
		try {
			// start thread service
			pool = Executors.newFixedThreadPool(threads);
			
			// prep tasks
			final SAMFileReader reader = new SAMFileReader(new File(bamFilename));
			final SAMSequenceDictionary dict = reader.getFileHeader().getSequenceDictionary();
			ArrayList< Callable<SAMRecord> > tasks = new ArrayList<Callable<SAMRecord>>();
			for (SAMSequenceRecord sq : dict.getSequences()) {
				tasks.add( new FindSeqByName(bamFilename, sq.getSequenceName(), sq.getSequenceLength(), readName));
			}
			
			reader.close();

			pool.invokeAll(tasks);
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			pool.shutdown();
		}
	}
	
	public FindSeqByName(String bamFilename, String chrom, int len, String readName) {
		readId = readName;
		this.chrom = chrom;
		chromLen = len;
		filename = bamFilename;
	}
	
	public SAMRecord call() {
		SAMRecord match = null;
		try {
			SAMFileReader sam = new SAMFileReader(new File(filename));
			SAMRecordIterator iter = sam.queryOverlapping(chrom, 1, chromLen);
	
			while (iter.hasNext()) {
				SAMRecord aln = iter.next();
				if (aln.getReadName().startsWith(readId)) {
					System.out.println( aln.getSAMString() );
					match = aln;
				}
			}
			iter.close();
			sam.close();
			
			if (match==null) {
				System.err.println("NOT FOUND IN: " + chrom);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}

		return match;
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		if (args.length==2) findSeqByName(args[0], args[1]);
		if (args.length==3) findSeqByName(args[0], args[1], Integer.parseInt(args[2]));
	}

}
