package au.edu.wehi.idsv.alignment;

import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqRecord;

import java.io.IOException;
import java.util.ArrayDeque;
import java.util.Queue;

/**
 * Performs full Smith-Waterman align across an entire chromosome. Only suitable for very small chromosomes.
 *
 */
public class SmithWatermanStreamingAligner implements StreamingAligner {
	private final Aligner aligner;
	private final int referenceIndex;
	private final ReferenceLookup reference;
	private final Queue<SAMRecord> buffer = new ArrayDeque<>();
	private final SAMFileHeader header;

	public SmithWatermanStreamingAligner(Aligner aligner, ReferenceLookup reference, int referenceIndex) {
		this.aligner = aligner;
		this.referenceIndex = referenceIndex;
		this.reference = reference;
		this.header = new SAMFileHeader();
		this.header.setSequenceDictionary(this.reference.getSequenceDictionary());
	}

	@Override
	public void asyncAlign(FastqRecord fq) throws IOException {
		byte[] ref = reference.getSequence(reference.getSequenceDictionary().getSequence(referenceIndex).getSequenceName()).getBases();
		Alignment alignment = aligner.align_smith_waterman(fq.getReadBases(), ref);
		SAMRecord r = SAMRecordUtil.createSAMRecord(header, fq, false);
		r.setReferenceIndex(referenceIndex);
		r.setAlignmentStart(alignment.getStartPosition() + 1);
		r.setCigarString(alignment.getCigar());
		buffer.add(r);
	}

	@Override
	public void flush() throws IOException {

	}

	@Override
	public int processedAlignmentRecords() {
		return buffer.size();
	}

	@Override
	public int outstandingAlignmentRecord() {
		return 0;
	}

	@Override
	public SAMRecord getAlignment() {
		return buffer.remove();
	}

	@Override
	public void close() throws IOException {

	}
}
