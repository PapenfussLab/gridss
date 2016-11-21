package au.edu.wehi.idsv.alignment;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang.ArrayUtils;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.GenomicProcessingContext;
import au.edu.wehi.idsv.SplitReadIdentificationHelper;
import au.edu.wehi.idsv.sam.ChimericAlignment;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.SequenceUtil;

/**
 * Performs full Smith-Waterman align across an entire chromosome. Only suitable for very small chromosomes.
 *
 */
public class StubFastqAligner implements FastqAligner {
	private final Map<String, ChimericAlignment> map = new HashMap<>();
	private final Map<String, SAMRecord> nameLookup = new HashMap<>();
	private final GenomicProcessingContext context;
	public StubFastqAligner(GenomicProcessingContext context) {
		this.context = context;
	}
	public StubFastqAligner align(SAMRecord r, BreakendDirection direction, int referenceIndex, int pos, boolean isNegativeStrand, String cigar) {
		List<FastqRecord> srs = SplitReadIdentificationHelper.getSplitReadRealignments(r, false);
		FastqRecord fqr = srs.get(0);
		if (srs.size() == 2) {
			if (direction == BreakendDirection.Forward ^ r.getReadNegativeStrandFlag()) {
				srs.get(1);
			}
		}
		map.put(fqr.getReadHeader(), new ChimericAlignment(r.getHeader().getSequenceDictionary().getSequence(referenceIndex).getSequenceName(),
				pos, isNegativeStrand, TextCigarCodec.decode(cigar), 40, 0));
		nameLookup.put(fqr.getReadHeader(), r);
		return this;
	}
	public StubFastqAligner align(SAMRecord r, int referenceIndex, int pos, boolean isNegativeStrand, String cigar) {
		List<FastqRecord> srs = SplitReadIdentificationHelper.getSplitReadRealignments(r, false);
		if(srs.size() != 1) throw new IllegalArgumentException("Need direction if source record contains clips on both sides");
		align(r, SAMRecordUtil.getStartClipLength(r) > 0 ? BreakendDirection.Backward : BreakendDirection.Forward, referenceIndex, pos, isNegativeStrand, cigar);
		return this;
	}
	@Override
	public void align(File fastq, File output, File reference, int threads) throws IOException {
		SAMFileHeader header = context.getBasicSamHeader();
		try (FastqReader reader = new FastqReader(fastq)) {
			try (SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, output)) {
				for (FastqRecord fqr : reader) {
					SAMRecord source = nameLookup.get(fqr.getReadString());
					SAMRecord r = new SAMRecord(header);
					r.setReadName(fqr.getReadHeader());
					r.setReadBases(fqr.getReadString().getBytes());
					r.setBaseQualities(SAMUtils.fastqToPhred(fqr.getBaseQualityString()));
					if (source == null) {
						r.setReadUnmappedFlag(true);
					} else {
						ChimericAlignment aln = map.get(source);
						r.setReferenceName(aln.rname);
						r.setAlignmentStart(aln.pos);
						r.setCigar(aln.cigar);
						if (aln.isNegativeStrand) {
							SequenceUtil.reverseComplement(r.getReadBases());
							ArrayUtils.reverse(r.getBaseQualities());
						}
					}
					writer.addAlignment(r);
				}
			}
		}
	}
}
