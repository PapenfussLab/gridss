package gridss;

import au.edu.wehi.idsv.ProgressLoggingSAMRecordIterator;
import au.edu.wehi.idsv.sam.ChimericAlignment;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.util.AsyncBufferedIterator;
import au.edu.wehi.idsv.util.ParallelTransformIterator;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import com.google.common.collect.Range;
import htsjdk.samtools.*;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.List;

@CommandLineProgramProperties(
		summary = "Exports unmapped sequences to fastq. " +
				"Only primary read alignments are exported.",
        oneLineSummary = "Exports unmapped sequences to fastq.",
        programGroup = gridss.cmdline.programgroups.DataConversion.class
)
public class UnmappedSequencesToFastq extends CommandLineProgram {
	private static final Log log = Log.getInstance(UnmappedSequencesToFastq.class);
	@Argument(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input SAM/BAM files")
	public List<File> INPUT;
	@Argument(shortName= StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output fastq file")
	public File OUTPUT;
	@Argument(doc="Minimum length of sequence export. " +
			"Generally speaking, very short soft clips are uninformative.", optional=true)
	public int MIN_SEQUENCE_LENGTH = 20;
	@Argument(doc="Include soft clipped bases. " +
			"For split read alignments, the largest contiguous sequence not aligned to the reference is used.", optional=true)
	public boolean INCLUDE_SOFT_CLIPPED_BASES = true;
	@Argument(doc="Ensure exported names are unique by suffixing with '/1' or '/2'", optional=true)
	public boolean UNIQUE_NAME = false;

	public static void main(String[] argv) {
		System.exit(new UnmappedSequencesToFastq().instanceMain(argv));
	}

	@Override
	protected int doWork() {
		for (File input : INPUT) {
			IOUtil.assertFileIsReadable(input);
		}
		IOUtil.assertFileIsWritable(OUTPUT);
		try (FastqWriter fqw = new FastqWriterFactory().newWriter(OUTPUT)) {
			for (File input : INPUT) {
				try (SamReader reader = SamReaderFactory.makeDefault().open(input)) {
					try (SAMRecordIterator rawIt = reader.iterator()) {
						ProgressLoggingSAMRecordIterator loggedIt = new ProgressLoggingSAMRecordIterator(rawIt, new ProgressLogger(log, 10000000));
						try (AsyncBufferedIterator<FastqRecord> it = new AsyncBufferedIterator(
								Iterators.filter(
										Iterators.transform(
												loggedIt,
												r -> getUnmappedFastqRecord(r)),
										r -> r != null),
								"toFastq")) {
							while (it.hasNext()) {
								FastqRecord fq = it.next();
								if (fq != null && fq.getReadLength() >= MIN_SEQUENCE_LENGTH) {
									fqw.write(fq);
								}
							}
						}
					}
				}
			}
		} catch (IOException e) {
			log.error(e);
			return -1;
		}
		return 0;
	}
	private FastqRecord getUnmappedFastqRecord(SAMRecord record) {
		FastqRecord fq = null;
		if (!record.isSecondaryOrSupplementary()) {
			byte[] bases = null;
			byte[] quals = null;
			if (record.getReadUnmappedFlag()) {
				bases = record.getReadBases();
				quals = record.getBaseQualities();
			} else if (INCLUDE_SOFT_CLIPPED_BASES && (SAMRecordUtil.getStartClipLength(record) >= MIN_SEQUENCE_LENGTH || SAMRecordUtil.getEndClipLength(record) >= MIN_SEQUENCE_LENGTH)) {
				// grab the largest contiguous subread not aligned
				List<ChimericAlignment> ca = Lists.newArrayList(new ChimericAlignment(record));
				ca.addAll(ChimericAlignment.getChimericAlignments(record));
				Range<Integer> mostUnaligned = Range.closed(0, 0);
				for (Range<Integer> r : ChimericAlignment.getUnalignedIntervals(ca).asRanges()) {
					if (r.upperEndpoint() - r.lowerEndpoint() > mostUnaligned.upperEndpoint() - mostUnaligned.lowerEndpoint()) {
						mostUnaligned = r;
					}
				}
				if (mostUnaligned.upperEndpoint() != mostUnaligned.lowerEndpoint()) {
					SAMRecordUtil.hardClipToN(record);
					bases = record.getReadBases().clone();
					quals = record.getReadBases().clone();
					if (record.getReadNegativeStrandFlag()) {
						SequenceUtil.reverseComplement(bases);
						SequenceUtil.reverseQualities(quals);
					}
					bases = Arrays.copyOfRange(bases, mostUnaligned.lowerEndpoint(), mostUnaligned.upperEndpoint());
					quals = Arrays.copyOfRange(quals, mostUnaligned.lowerEndpoint(), mostUnaligned.upperEndpoint());
				}
			}
			if (bases != null) {
				String uniqueName = record.getReadName();
				if (UNIQUE_NAME) {
					if (record.getFirstOfPairFlag()) {
						uniqueName += "/1";
					} else if (record.getSecondOfPairFlag()) {
						uniqueName += "/2";
					}
				}
				fq = new FastqRecord(
						uniqueName,
						new String(bases, StandardCharsets.UTF_8),
						null,
						SAMUtils.phredToFastq(quals));
			}
		}
		return fq;
	}
}
