package gridss;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;

import com.google.common.collect.Lists;

import au.edu.wehi.idsv.ReadPairConcordanceCalculator;
import au.edu.wehi.idsv.ReadPairConcordanceMethod;
import au.edu.wehi.idsv.metrics.InsertSizeDistribution;
import au.edu.wehi.idsv.util.AsyncBufferedIterator;
import gridss.filter.ClippedReadFilter;
import gridss.filter.ReadPairConcordanceFilter;
import gridss.filter.IndelReadFilter;
import gridss.filter.OneEndAnchoredReadFilter;
import gridss.filter.SplitReadFilter;
import gridss.filter.UnionAggregateFilter;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.filter.AlignedFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(
        usage = "Uses an external aligner to identify split reads by iterative alignment of soft clipped bases.",
        usageShort = "Converts soft clipped reads to split reads"
)
public class SoftClipsToSplitReads extends CommandLineProgram {
	private static final Log log = Log.getInstance(SoftClipsToSplitReads.class);
	private static final int ASYNC_BUFFERS = 2;
	private static final int ASYNC_BUFFER_SIZE = 300;
    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input file", optional=false)
    public File INPUT;
    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output file", optional=false)
    public File OUTPUT;
    @Option(doc="Minimum bases clipped. Generally, short read aligners are not able to align sequences shorter than 18-20 bases.", optional=true)
    public int MIN_CLIP_LENGTH = 15;
    @Override
	protected int doWork() {
		log.debug("Setting language-neutral locale");
    	java.util.Locale.setDefault(Locale.ROOT);
    	validateParameters();
    	
    	SamReaderFactory readerFactory = SamReaderFactory.make();
    	
    	List<File> aligned = Lists.newArrayList(INPUT);
    	int recordsWritten = 1;
    	while (recordsWritten > 0) {
    		// foreach SAMRecord r in current BAM
    			//foreach realignFastaRecords(r)
    				// write fastq
    				// recordsWritten++ 
    		// align fastq
    	}
    	// merge back into a SA record
    	SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
    	try (SamReader reader = readerFactory.open(INPUT)) {
			SAMFileHeader header = reader.getFileHeader();
			try (SAMRecordIterator it = reader.iterator()) {
				try (SAMFileWriter writer = writerFactory.makeSAMOrBAMWriter(header, true, OUTPUT)) {
					while (it.hasNext()) {
						SAMRecord r = it.next();
						List<SAMRecord> alignments = Lists.newArrayListWithExpectedSize(4);
						alignments.add(r);
						int expected = realignFastaRecords(r).size();
						for (SAMRecordIterator realignIt) {
							// find all alignments - they should be the first records
							while (isSameRead(r, realignIt.peek())) {
								SAMRecord realignedsc = realignIt.next();
								alignments.add(realignedsc);
							}
						}
						createSupplimentaryAlignments(alignments);
						writeBAM(r);
						alignments.remove(r);
						writeSortedBAM(supplimentary)
					}
				}
			}
		}
    	// if sort order == coordinate
    	// sort supplimentary bam file
    	// merge BAM and supplimentary BAM
    	try {
		} catch (IOException e) {
			log.error(e);
			return -1;
		}
    	return 0;
	}
    
	private void validateParameters() {
    	IOUtil.assertFileIsReadable(INPUT);
    	IOUtil.assertFileIsWritable(OUTPUT);
	}
	public static void main(String[] argv) {
        System.exit(new SoftClipsToSplitReads().instanceMain(argv));
    }
}
