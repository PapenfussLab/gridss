package au.edu.wehi.validation;

import au.edu.wehi.idsv.sam.ChimericAlignment;
import au.edu.wehi.idsv.sam.CigarUtil;
import htsjdk.samtools.*;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.tribble.bed.SimpleBEDFeature;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Validations deletion calls with long read support
 * @author Daniel Cameron
 *
 */
@CommandLineProgramProperties(
		summary = "Converts split read alignments to bed intervals",
		oneLineSummary = "Converts split read alignments to bed intervals",
		programGroup = gridss.cmdline.programgroups.DataConversion.class)
public class BamToBed extends CommandLineProgram {
    @Argument(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, optional=false)
    public File INPUT;
    @Argument(doc="Split read support", shortName="SR", optional=true)
    public File SPLIT_READS;
    @Argument(doc="Indel support for reads spanning", shortName="SP", optional=true)
    public File SPANNING_READS;
    @Argument(doc="Only report deletions", shortName="OD", optional=true)
    public boolean ONLY_DELETIONS = true;
    @Argument(doc="Minimum percentage of spplit read bases aligned to each location. Chimeric alignments aligning few than this portion of the read will be ignored.", shortName="MAP", optional=true)
    public double MINIMUM_ALIGNED_PORTION = 0.25;
    @Argument(doc="Minimum mapq to be considered a valid alignment.", shortName="MQ", optional=true)
    public int MAPQ = 1;
    @Argument(doc="Ignore indels smaller than this size", shortName="MIS", optional=true)
    public int MIN_FINAL_INDEL_SIZE = 50;
    @Argument(doc="Merge indels with few than this many reads aligned. This is useful for PacBio reads alignments that contain fragmented indel calls.", shortName="MIMB", optional=true)
    public int MAX_INTERVENING_MAPPED_BASES = 0;
    @Argument(doc="Minimum size of partial indel. This is useful for PacBio reads alignments that contain fragmented indel calls.", shortName="MCIS", optional=true)
    public int MIN_COMPONENT_INDEL_SIZE = 5;
    protected int doWork() {
    	if (!ONLY_DELETIONS) {
    		System.err.append("Processing non-deletion events not yet implemented");
    		return 1;
    	}
        try {
        	extract(INPUT, SPLIT_READS, SPANNING_READS, MINIMUM_ALIGNED_PORTION, MAPQ, MIN_FINAL_INDEL_SIZE, MAX_INTERVENING_MAPPED_BASES, MIN_COMPONENT_INDEL_SIZE);
        } catch (Exception e) {
			e.printStackTrace();
			return 1;
		}
        return 0;
    }
    public static void extract(File input, File splitReads, File spanningReads, final double minLengthPortion, final int minMapq, final int minIndelSize, final int maxMappedBases, final int minIndelComponentSize) throws IOException {
    	SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT);
		SamReader reader = factory.open(input);
		
    	BufferedWriter srwriter = null;
    	if (splitReads != null) {
    		srwriter = new BufferedWriter(new FileWriter(splitReads));
    	}
    	BufferedWriter spwriter = null;
    	if (spanningReads != null) {
    		spwriter = new BufferedWriter(new FileWriter(spanningReads));
    	}
		SAMRecordIterator it = reader.iterator();
		while (it.hasNext()) {
			final SAMRecord r = it.next();
			if (r.isSecondaryOrSupplementary()) continue;
			if (splitReads != null) {
				for (SimpleBEDFeature f : ChimericAlignment.getChimericAlignments(r).stream().map(ca -> getSplitReadDeletion(r, ca, minLengthPortion, minMapq, minIndelSize))
						.filter(f -> f != null)
						.collect(Collectors.toList())) {
					writeBed(srwriter, f);
				}
			}
			if (spanningReads != null) {
				for (SimpleBEDFeature f : getSpanningDeletion(r, minMapq, minIndelSize, maxMappedBases, minIndelComponentSize)) {
					writeBed(spwriter, f);
				}
			}
		}
		CloserUtil.close(srwriter);
		CloserUtil.close(spwriter);
    }
    private static void writeBed(Writer writer, SimpleBEDFeature feature) throws IOException {
    	writer.write(feature.getContig());
		writer.write('\t');
		writer.write(Integer.toString(feature.getStart() + 1));
		writer.write('\t');
		writer.write(Integer.toString(feature.getEnd() + 1));
		writer.write('\n');
    }
    public static List<SimpleBEDFeature> getSpanningDeletion(final SAMRecord r, final int minMapq, final int minIndelSize, final int maxMappedBases, final int minIndelComponentSize) {
    	List<SimpleBEDFeature> list = new ArrayList<SimpleBEDFeature>();
    	if (r.getMappingQuality() < minMapq) return list;
    	List<CigarElement> ce = r.getCigar().getCigarElements();
    	int offset = 0;
    	for (int i = 0; i < ce.size(); i++) {
    		CigarElement e = ce.get(i);
    		if (e.getOperator() == CigarOperator.D) {
    			if (maxMappedBases > 0) {
    				throw new RuntimeException("Indel merging via MAX_INTERVENING_MAPPED_BASES not yet implemented.");
    			}
    			if (e.getLength() > minIndelSize) {
    				list.add(new SimpleBEDFeature(r.getAlignmentStart() + offset, r.getAlignmentStart() + offset + e.getLength(), r.getReferenceName()));
    			}
    		}
    		if (e.getOperator().consumesReferenceBases()) {
    			offset += e.getLength();
    		}
    	}
    	return list;
    }
    public static SimpleBEDFeature getSplitReadDeletion(SAMRecord r, ChimericAlignment ca, double minLengthPortion, int minMapq, int minIndelSize) {
    	if (r.getMappingQuality() < minMapq) return null;
    	if (ca.mapq < minMapq) return null;
    	if (!ca.rname.equals(r.getReferenceName())) return null;
    	if (ca.isNegativeStrand != r.getReadNegativeStrandFlag()) return null;
    	if (r.getCigar().getReferenceLength() < r.getReadLength() * minLengthPortion) return null;
    	if (ca.cigar.getReferenceLength() < r.getReadLength() * minLengthPortion) return null;
    	int caStartClipLength = CigarUtil.getStartClipLength(ca.cigar.getCigarElements());
    	int rStartClipLength = CigarUtil.getStartClipLength(r.getCigar().getCigarElements());
    	if (r.getAlignmentEnd() < ca.pos && caStartClipLength > rStartClipLength) {
    		if (ca.pos - r.getAlignmentEnd() >= minIndelSize) { 
    			return new SimpleBEDFeature(r.getAlignmentEnd(), ca.pos, r.getReferenceName());
    		}
    	}
    	int endpos = ca.pos + ca.cigar.getReferenceLength() - 1;
    	if (endpos < r.getAlignmentStart() && caStartClipLength < rStartClipLength) {
    		if (r.getAlignmentStart() - endpos >= minIndelSize) {
    			return new SimpleBEDFeature(endpos, r.getAlignmentStart(), r.getReferenceName());
    		}
    	}
    	return null;
    }
	public static void main(String[] argv) {
	    System.exit(new BamToBed().instanceMain(argv));
	}
}
