package au.edu.wehi.idsv.sam;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;

import java.util.Iterator;

/**
 * Cleans CIGARs so they comply with the SAM specifications
 * @author cameron.d
 *
 */
public class SAMRecordCigarCleaningIterator implements CloseableIterator<SAMRecord> {
	private static final Log log = Log.getInstance(SAMRecordCigarCleaningIterator.class);
	private final Iterator<SAMRecord> underlying;
	public SAMRecordCigarCleaningIterator(Iterator<SAMRecord> underlying) {
		this.underlying = underlying;
	}
	@Override
	public boolean hasNext() {
		return underlying.hasNext();
	}

	@Override
	public SAMRecord next() {
		return cleanCigar(underlying.next());
	}

	private SAMRecord cleanCigar(SAMRecord next) {
		Cigar cigar = next.getCigar();
		if (cigar != null && cigar.getCigarElements().size() > 0) {
			Cigar newCigar = new Cigar(CigarUtil.clean(cigar.getCigarElements(), true));
			if (!cigar.equals(newCigar)) {
				log.warn("Cigar %s of read %s is not valid.", cigar, next.getReadName());
				next.setCigar(newCigar);
			}
		}
		return next;
	}
	@Override
	public void close() {
		CloserUtil.close(underlying);
	}

}
