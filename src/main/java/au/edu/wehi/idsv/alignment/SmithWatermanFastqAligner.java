package au.edu.wehi.idsv.alignment;

import htsjdk.samtools.*;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFile;

import java.io.File;
import java.io.IOException;

/**
 * Performs full Smith-Waterman align across an entire chromosome. Only suitable for very small chromosomes.
 *
 */
public class SmithWatermanFastqAligner implements FastqAligner {
	private final Aligner aligner;
	private final int referenceIndex;
	public SmithWatermanFastqAligner(Aligner aligner, int referenceIndex) {
		this.aligner = aligner;
		this.referenceIndex = referenceIndex;
	}
	@Override
	public void align(File fastq, File output, File reference, int threads) throws IOException {
		try (ReferenceSequenceFile ref = new IndexedFastaSequenceFile(reference)) {
			SAMFileHeader header = new SAMFileHeader();
			header.setSequenceDictionary(ref.getSequenceDictionary());
			byte[] bases = ref.getSequence(ref.getSequenceDictionary().getSequence(referenceIndex).getSequenceName()).getBases();
			try (FastqReader reader = new FastqReader(fastq)) {
				try (SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, output)) {
					for (FastqRecord fqr : reader) {
						Alignment aln = aligner.align_smith_waterman(fqr.getReadString().getBytes(), bases);
						SAMRecord r = new SAMRecord(header);
						r.setReadName(fqr.getReadName());
						r.setReferenceIndex(referenceIndex);
						r.setAlignmentStart(aln.getStartPosition() + 1);
						r.setCigarString(aln.getCigar());
						r.setReadBases(fqr.getReadString().getBytes());
						r.setBaseQualities(SAMUtils.fastqToPhred(fqr.getBaseQualityString()));
						writer.addAlignment(r);
					}
				}
			}
		}
	}
}
