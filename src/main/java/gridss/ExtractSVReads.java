package gridss;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import au.edu.wehi.idsv.FileSystemContext;
import au.edu.wehi.idsv.ReadPairConcordanceCalculator;
import au.edu.wehi.idsv.ReadPairConcordanceMethod;
import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.sam.ChimericAlignment;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.util.FileHelper;
import gridss.analysis.CollectStructuralVariantReadMetrics;
import gridss.cmdline.ProcessStructuralVariantReadsCommandLineProgram;
import gridss.filter.ClippedReadFilter;
import gridss.filter.FixedFilter;
import gridss.filter.IndelReadFilter;
import gridss.filter.OneEndAnchoredReadFilter;
import gridss.filter.ReadPairConcordanceFilter;
import gridss.filter.SplitReadFilter;
import gridss.filter.UnionAggregateFilter;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.filter.AlignedFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.Log;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;

@CommandLineProgramProperties(
        usage = "Extracts reads and read pairs supporting putative structural variations. "
        		+ "If the input file is queryname sorted, a multi-mapping aware extraction is performed "
        		+ "and reads/read pairs are only extracted when all alignments are consistent with the "
        		+ "presence of of a structural variant.",
        usageShort = "Extracts reads and read pairs supporting putative structural variations."
)
public class ExtractSVReads extends ProcessStructuralVariantReadsCommandLineProgram {
	private static final Log log = Log.getInstance(ExtractSVReads.class);
    @Option(shortName="MO", doc="Output file containing SV metrics", optional=true)
    public File METRICS_OUTPUT;
    private CollectStructuralVariantReadMetrics metricsCollector;
    private File tmpoutput;
    private SAMFileWriter writer;
    private SamRecordFilter readfilter;
    private SamRecordFilter pairfilter;
    private int count;
    @Override
    protected void setup(SAMFileHeader header, File samFile) {
    	if (METRICS_OUTPUT != null) {
    		metricsCollector = new CollectStructuralVariantReadMetrics();
    		copyInput(metricsCollector);
    		metricsCollector.OUTPUT = METRICS_OUTPUT;
    		metricsCollector.setup(header, samFile);
    	}
    	SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
    	if (header.getSortOrder() != SortOrder.queryname) {
			log.info("Not considering multiple read alignments as the input file is not queryname sorted.");
		}
    	tmpoutput = FileSystemContext.getWorkingFileFor(OUTPUT, "gridss.tmp.ExtractSVReads.");
    	writer = writerFactory.makeSAMOrBAMWriter(header, true, tmpoutput);
    	
    	IndelReadFilter indelFilter = new IndelReadFilter(INDELS ? MIN_INDEL_SIZE : Integer.MAX_VALUE);
		ClippedReadFilter softClipFilter = new ClippedReadFilter(CLIPPED ? MIN_CLIP_LENGTH : Integer.MAX_VALUE); 
		SplitReadFilter splitReadFilter = new SplitReadFilter();
		AlignedFilter unmappedFilter = new AlignedFilter(false);
		OneEndAnchoredReadFilter oeaFilter = new OneEndAnchoredReadFilter();
		ReadPairConcordanceFilter dpFilter = new ReadPairConcordanceFilter(getReadPairConcordanceCalculator(), false, true);
		List<SamRecordFilter> readfilters = new ArrayList<>();
		readfilters.add(indelFilter);
		readfilters.add(softClipFilter);
		if (SPLIT) readfilters.add(splitReadFilter);
		if (UNMAPPED_READS) readfilters.add(unmappedFilter);
		readfilter = new UnionAggregateFilter(readfilters);
		
		List<SamRecordFilter> pairfilters = new ArrayList<>();
		if (SINGLE_MAPPED_PAIRED) pairfilters.add(oeaFilter);
		if (DISCORDANT_READ_PAIRS) pairfilters.add(dpFilter);
		pairfilter = new UnionAggregateFilter(pairfilters);
		if (!SINGLE_MAPPED_PAIRED && !DISCORDANT_READ_PAIRS) {
			pairfilter = new FixedFilter(true);
		}
		count = 0;
    }
	public static boolean[] hasReadAlignmentConsistentWithReference(List<SAMRecord> records) {
		boolean[] consistent = new boolean[2];
		for (SAMRecord r : records) {
			int segmentIndex = SAMRecordUtil.getSegmentIndex(r);
			if (consistent.length <= segmentIndex) {
				consistent = Arrays.copyOf(consistent, segmentIndex + 1);
			}
			if (consistent[segmentIndex]) continue; // no need to check this segment any further
			if (isFullyMapped(r)) {
				consistent[segmentIndex] = true;
			}
			if (hasFullyMappedSplit(r)) {
				consistent[segmentIndex] = true;
			}
		}
		return consistent;
	}
	private static boolean hasFullyMappedSplit(SAMRecord r) {
		for (ChimericAlignment ca : ChimericAlignment.getChimericAlignments(r)) {
			if (isFullyMapped(ca.cigar.getCigarElements())) {
				return true;
			}
		}
		return false;
	}
	private static boolean isFullyMapped(SAMRecord r) {
		if (r.getReadUnmappedFlag()) return false;
		Cigar cigar = r.getCigar();
		if (cigar == null) return false;
		return isFullyMapped(cigar.getCigarElements());
	}
	private static boolean isFullyMapped(List<CigarElement> cigar) {
		for (CigarElement ce : cigar) {
			switch (ce.getOperator()) {
				case M:
				case EQ:
				case X:
				case P:
					break;
				default:
					return false;
			}
		}
		return true;
	}
	public static boolean hasReadPairingConsistentWithReference(ReadPairConcordanceCalculator rpcc, List<SAMRecord> records) {
		if (rpcc == null) return false;
		for (SAMRecord r1 : records) {
			if (mateIsConsistentWithReference(rpcc, r1)) {
				return true;
			}
			if (r1.getReadUnmappedFlag() || !r1.getReadPairedFlag() || !r1.getFirstOfPairFlag()) continue;
			for (SAMRecord r2 : records) {
				if (r2.getReadUnmappedFlag() || !r2.getReadPairedFlag() || !r2.getSecondOfPairFlag()) continue;
				if (rpcc.isConcordant(r1, r2)) {
					return true;
				}
			}
		}
		return false;
	}
	public static boolean mateIsConsistentWithReference(ReadPairConcordanceCalculator rpcc, SAMRecord record) {
		if (rpcc == null) return false;
		if (record.getReadUnmappedFlag()) return false;
		if (!record.getReadPairedFlag()) return false;
		if (record.getMateUnmappedFlag()) return false;
		if (record.getSupplementaryAlignmentFlag()) {
			// Use the primary alignment for supp alignment
			record = SAMRecordUtil.clone(record);
			List<ChimericAlignment> calist = ChimericAlignment.getChimericAlignments(record);
			if (calist.size() > 0) {
				ChimericAlignment ca = calist.get(0);
				record.setReferenceName(ca.rname);
				record.setAlignmentStart(ca.pos);
				record.setReadNegativeStrandFlag(ca.isNegativeStrand);
				record.setCigar(ca.cigar);
				record.setMappingQuality(ca.mapq);
			}
		}
		return rpcc.isConcordant(record, null);
	}
	@Override
	protected String[] customCommandLineValidation() {
		if (READ_PAIR_CONCORDANCE_METHOD == ReadPairConcordanceMethod.PERCENTAGE && INSERT_SIZE_METRICS == null) {
			return new String[] { "INSERT_SIZE_METRICS is required when using percentage based read pair concordance" };
		}
		return super.customCommandLineValidation();
	}
	public static void main(String[] argv) {
        System.exit(new ExtractSVReads().instanceMain(argv));
    }
	@Override
	protected void acceptFragment(List<SAMRecord> records, ReferenceLookup lookup) {
		boolean hasConsistentReadPair = hasReadPairingConsistentWithReference(getReadPairConcordanceCalculator(), records);
		boolean[] hasConsistentReadAlignment = hasReadAlignmentConsistentWithReference(records);
		if (metricsCollector != null) {
			metricsCollector.acceptFragment(records, lookup);
		}
		for (SAMRecord r : records) {
			if ((!readfilter.filterOut(r) && !hasConsistentReadAlignment[SAMRecordUtil.getSegmentIndex(r)])
					|| (!pairfilter.filterOut(r) && !hasConsistentReadPair)) {
				writer.addAlignment(r);
				count++;
			} else {
				// ignore remaining reads
			}
		}
	}
	@Override
	protected void finish() {
		writer.close();
		try {
			FileHelper.move(tmpoutput, OUTPUT, true);
			log.info(String.format("Extracted %d reads from %s", count, INPUT));
		} catch (IOException e) {
			log.error(e);
			throw new RuntimeException(e);
		}
		if (METRICS_OUTPUT != null) {
			metricsCollector.finish();
		}
	}
}
