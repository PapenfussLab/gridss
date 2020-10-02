package gridss;

import au.edu.wehi.idsv.FileSystemContext;
import au.edu.wehi.idsv.ReadPairConcordanceCalculator;
import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.sam.ChimericAlignment;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.util.FileHelper;
import gridss.cmdline.ProcessStructuralVariantReadsCommandLineProgram;
import gridss.filter.*;
import htsjdk.samtools.*;
import htsjdk.samtools.filter.AlignedFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

@CommandLineProgramProperties(
		summary = "Extracts reads and read pairs supporting putative structural variations. "
        		+ "If the input file is queryname sorted, a multi-mapping aware extraction is performed "
        		+ "and reads/read pairs are only extracted when all alignments are consistent with the "
        		+ "presence of of a structural variant.",
        oneLineSummary = "Extracts reads and read pairs supporting putative structural variations.",
        programGroup = picard.cmdline.programgroups.ReadDataManipulationProgramGroup.class
)
public class ExtractSVReads extends ProcessStructuralVariantReadsCommandLineProgram {
	private static final Log log = Log.getInstance(ExtractSVReads.class);
    private File tmpoutput;
    private SAMFileWriter writer;
    private SamRecordFilter readfilter;
    private SamRecordFilter pairfilter;
    private int count;
    @Override
    protected void setup(SAMFileHeader header, File samFile) {
    	SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
    	tmpoutput = gridss.Defaults.OUTPUT_TO_TEMP_FILE ? FileSystemContext.getWorkingFileFor(OUTPUT, "gridss.tmp.ExtractSVReads.") : OUTPUT;
    	writer = writerFactory.makeSAMOrBAMWriter(header, true, tmpoutput);
    	
    	IndelReadFilter indelFilter = new IndelReadFilter(INDELS ? MIN_INDEL_SIZE : Integer.MAX_VALUE);
		ClippedReadFilter softClipFilter = new ClippedReadFilter(CLIPPED ? MIN_CLIP_LENGTH : Integer.MAX_VALUE); 
		SplitReadFilter splitReadFilter = new SplitReadFilter();
		AlignedFilter unmappedFilter = new AlignedFilter(false);
		OneEndAnchoredReadFilter oeaFilter = new OneEndAnchoredReadFilter();
		ReadPairConcordanceFilter dpFilter = getReadPairConcordanceCalculator() != null ? new ReadPairConcordanceFilter(getReadPairConcordanceCalculator(), false, true) : null;
		List<SamRecordFilter> readfilters = new ArrayList<>();
		readfilters.add(indelFilter);
		readfilters.add(softClipFilter);
		if (SPLIT) readfilters.add(splitReadFilter);
		if (UNMAPPED_READS) readfilters.add(unmappedFilter);
		readfilter = new UnionAggregateFilter(readfilters);
		
		List<SamRecordFilter> pairfilters = new ArrayList<>();
		if (SINGLE_MAPPED_PAIRED) pairfilters.add(oeaFilter);
		if (dpFilter != null && DISCORDANT_READ_PAIRS) pairfilters.add(dpFilter);
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
		return rpcc.isConcordant(primaryAlignmentForSupplementary(record), null);
	}
	private static SAMRecord primaryAlignmentForSupplementary(SAMRecord r) {
		if (r.getSupplementaryAlignmentFlag()) {
			SAMRecord record = SAMRecordUtil.clone(r);
			List<ChimericAlignment> calist = ChimericAlignment.getChimericAlignments(r);
			if (calist.size() > 0) {
				ChimericAlignment ca = calist.get(0);
				record.setReferenceName(ca.rname);
				record.setAlignmentStart(ca.pos);
				record.setReadNegativeStrandFlag(ca.isNegativeStrand);
				record.setCigar(ca.cigar);
				record.setMappingQuality(ca.mapq);
				record.setReadBases(SAMRecord.NULL_SEQUENCE);
				record.setBaseQualities(SAMRecord.NULL_QUALS);
				return record;
			}
		}
		return r;
	}
	public static void main(String[] argv) {
        System.exit(new ExtractSVReads().instanceMain(argv));
    }
	public boolean[] shouldExtract(List<SAMRecord> records, ReferenceLookup lookup) {
		boolean hasConsistentReadPair = hasReadPairingConsistentWithReference(getReadPairConcordanceCalculator(), records);
		boolean[] hasConsistentReadAlignment = hasReadAlignmentConsistentWithReference(records);
		boolean[] extract = new boolean[records.size()];
		for (int i = 0; i < records.size(); i++) {
			SAMRecord r = records.get(i);
			extract[i] = !hasConsistentReadAlignment[SAMRecordUtil.getSegmentIndex(r)] && !readfilter.filterOut(r);
			// supp records should use the primary alignment when considering concordance
			extract[i] |= !hasConsistentReadPair && !pairfilter.filterOut(primaryAlignmentForSupplementary(r));
			extract[i] &= (!r.getDuplicateReadFlag() || INCLUDE_DUPLICATES);
		}
		return extract;
	}
	@Override
	protected void acceptFragment(List<SAMRecord> records, ReferenceLookup lookup) {
		boolean[] extract = shouldExtract(records, lookup);
		for (int i = 0; i < records.size(); i++) {
			SAMRecord r = records.get(i);
			if (extract[i]) {
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
			if (tmpoutput != OUTPUT) {
				FileHelper.move(tmpoutput, OUTPUT, true);
			}
			log.info(String.format("Extracted %d reads from %s", count, INPUT));
		} catch (IOException e) {
			log.error(e);
			throw new RuntimeException(e);
		}
	}
}
