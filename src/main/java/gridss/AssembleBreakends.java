package gridss;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;

import org.apache.commons.lang.NotImplementedException;

import com.google.common.collect.Iterators;

import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.ReadPairConcordanceCalculator;
import au.edu.wehi.idsv.ReadPairConcordanceMethod;
import au.edu.wehi.idsv.SAMRecordAssemblyEvidence;
import au.edu.wehi.idsv.debruijn.positional.DirectedPositionalAssembler;
import au.edu.wehi.idsv.debruijn.positional.PositionalAssembler;
import au.edu.wehi.idsv.util.AsyncBufferedIterator;
import gridss.analysis.InsertSizeDistribution;
import gridss.filter.ClippedReadFilter;
import gridss.filter.IndelReadFilter;
import gridss.filter.OneEndAnchoredReadFilter;
import gridss.filter.ReadPairConcordanceFilter;
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
import htsjdk.variant.variantcontext.writer.BCF2FieldEncoder.AtomicInt;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(
        usage = "Extracts reads and read pairs supporting putative structural variations.",
        usageShort = "Extracts reads and read pairs supporting putative structural variations."
)
public class AssembleBreakends extends CommandLineProgram {
	private static final Log log = Log.getInstance(AssembleBreakends.class);
    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input file", optional=false)
    public List<File> INPUT;
    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output file containing subset of input", optional=false)
    public File OUTPUT;
    @Override
	protected int doWork() {
		log.debug("Setting language-neutral locale");
    	java.util.Locale.setDefault(Locale.ROOT);
    	validateParameters();
    	ProcessingContext pc = null;
    	SAMFileHeader header = pc.getBasicSamHeader();
    }
    
	public static void extract(
			final Iterator<SAMRecord> rawit,
			final SAMFileWriter writer,
			final int minIndelSize,
			final int minClipLength,
			final boolean includeSplitReads,
			final boolean includeOEA,
			final boolean includeDP,
			final ReadPairConcordanceCalculator rpcc,
			final boolean includeUnmapped) throws IOException {
		ProgressLogger progress = new ProgressLogger(log);
		List<SamRecordFilter> filters = new ArrayList<>();
		filters.add(new IndelReadFilter(minIndelSize));
		filters.add(new ClippedReadFilter(minClipLength));
		if (includeSplitReads) filters.add(new SplitReadFilter());
		if (includeOEA) filters.add(new OneEndAnchoredReadFilter());
		if (includeDP) filters.add(new ReadPairConcordanceFilter(rpcc, false, true));
		if (includeUnmapped) filters.add(new AlignedFilter(false));
		
		UnionAggregateFilter filter = new UnionAggregateFilter(filters);
		try (CloseableIterator<SAMRecord> it = new AsyncBufferedIterator<SAMRecord>(rawit, "raw")) {
			while (it.hasNext()) {
				SAMRecord r = it.next();
				progress.record(r);
				if (!filter.filterOut(r)) {
					writer.addAlignment(r);
				}
			}
		}
	}
	private void validateParameters() {
    	IOUtil.assertFileIsReadable(INPUT);
    	IOUtil.assertFileIsWritable(OUTPUT);
	}
	public static void main(String[] argv) {
        System.exit(new AssembleBreakends().instanceMain(argv));
    }
}
