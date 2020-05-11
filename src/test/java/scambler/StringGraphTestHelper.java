package scambler;

import au.edu.wehi.idsv.TestHelper;
import htsjdk.samtools.SAMRecord;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

public class StringGraphTestHelper extends TestHelper {
	public static List<SAMRecord> overlapping(byte[] seq, int readLength, int seqlength, int stride, boolean isCircular) {
		String seq1 = S(seq).substring(0, seqlength);
		String seq2 = S(seq).substring(seqlength, 2 * seqlength);
		if (isCircular) {
			seq2 = seq1;
		}
		List<SAMRecord> reads = new ArrayList<>();
		for (int i = 0; i < seq1.length(); i += stride) {
			String perfect_read = (seq1 + seq2).substring(i, i + readLength);
			SAMRecord r = new SAMRecord(getHeader());
			r.setReadName(Integer.toString(i));
			r.setReadBases(B(perfect_read));
			r.setReferenceIndex(0);
			r.setAlignmentStart(i);
			reads.add(r);
		}
		return reads;
	}
	public static List<SAMRecord> overlapping(int readCount, int readLength, int stride) {
		List<SAMRecord> reads = new ArrayList<>();
		for (int i = 0; i < readCount; i++) {
			SAMRecord r = new SAMRecord(getHeader());
			r.setReadName(new String(new byte[] {(byte)('A' + i)}));
			r.setReadBases(B(S(RANDOM).substring(i * stride, i * stride + readLength)));
			r.setCigarString(String.format("%dM", readLength));
			r.setReferenceIndex(3);
			r.setAlignmentStart(i * stride);
			reads.add(r);
		}
		return reads;
	}
	@Test
	public void test_overlapping() {
		List<SAMRecord> records = overlapping(8, 1, 1);
	}
}