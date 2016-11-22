package au.edu.wehi.idsv.util;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

import au.edu.wehi.idsv.ProgressLoggingSAMRecordIterator;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(
        usage = "Updates paired end reads without a mate to single-ended reads",  
        usageShort = "Updates paired end reads without a mate to single-ended reads"
)
public class FixMissingMateInformation extends picard.cmdline.CommandLineProgram {
	private Log log = Log.getInstance(FixMissingMateInformation.class);
	@Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="input file")
    public File INPUT;
	@Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="output file")
    public File OUTPUT;
	@Override
	protected int doWork() {
		SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
		SamReader input = factory.open(INPUT);
		Iterator<SAMRecord> it = new AsyncBufferedIterator<SAMRecord>(input.iterator(), 2, 128);
		SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(input.getFileHeader(), true, OUTPUT);
		it = new ProgressLoggingSAMRecordIterator(it, new ProgressLogger(log));
		PeekingIterator<SAMRecord> pit = Iterators.peekingIterator(it);
		while (pit.hasNext()) {
			List<SAMRecord> readlist = getRecordsForName(pit);
			fixPairingFlags(readlist);
			for (SAMRecord r : readlist) {
				out.addAlignment(r);
			}
		}
		out.close();
		return 0;
	}
	private void fixPairingFlags(List<SAMRecord> readlist) {
		boolean foundFirst = false;
		boolean foundSecond = false;
		for (SAMRecord r : readlist) {
			foundFirst |= r.getFirstOfPairFlag();
			foundSecond |= r.getSecondOfPairFlag();
		}
		if (!foundFirst || !foundSecond) {
			for (SAMRecord r : readlist) {
				r.setFirstOfPairFlag(false);
				r.setSecondOfPairFlag(false);
				r.setReadPairedFlag(false);
				r.setMateAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
				r.setMateNegativeStrandFlag(false);
				r.setMateReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
				r.setMateUnmappedFlag(false);
			}
		}
	}
	private List<SAMRecord> getRecordsForName(PeekingIterator<SAMRecord> it) {
		List<SAMRecord> list = new ArrayList<SAMRecord>(2);
		SAMRecord r = it.next();
		String name = r.getReadName();
		list.add(r);
		while (it.hasNext() && name != null && name.equals(it.peek().getReadName())) {
			list.add(it.next());
		}
		return list;
	}
	public static void main(String[] argv) {
        System.exit(new FixMissingMateInformation().instanceMain(argv));
    }
}
