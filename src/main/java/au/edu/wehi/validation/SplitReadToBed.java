package au.edu.wehi.validation;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.tribble.bed.SimpleBEDFeature;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import au.edu.wehi.idsv.sam.ChimericAlignment;
import au.edu.wehi.idsv.sam.CigarUtil;

/**
 * Validations deletion calls with long read support
 * @author Daniel Cameron
 *
 */
@CommandLineProgramProperties(
		usage = "Converts split read alignments to bed intervals",
		usageShort = "Converts split read alignments to bed intervals")
public class SplitReadToBed extends CommandLineProgram {
    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME)
    public File OUTPUT;
    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File INPUT;
    @Option(doc="Only report deletions", shortName="OD", optional=true)
    public boolean ONLY_DELETIONS = true;
    @Option(doc="Minimum percentage of read aligned to each location. Chimeric alignments aligning few than this portion of the read will be ignored.", shortName="MAP", optional=true)
    public double MINIMUM_ALIGNED_PORTION = 0.25;
    @Option(doc="Minimum mapq to be considered a valid alignment.", shortName="MQ", optional=true)
    public int MAPQ = 1;
    protected int doWork() {
    	if (!ONLY_DELETIONS) {
    		System.err.append("Processing non-deletion events not yet implemented");
    		return 1;
    	}
        try {
        	extractSplitReads(INPUT, OUTPUT, MINIMUM_ALIGNED_PORTION, MAPQ);
        } catch (Exception e) {
			e.printStackTrace();
			return 1;
		}
        return 0;
    }
    public static void extractSplitReads(File file, File output, double minLengthPortion, int minMapq) throws IOException {
    	SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT);
		SamReader reader = factory.open(file);
    	BufferedWriter writer = new BufferedWriter(new FileWriter(output));
		SAMRecordIterator it = reader.iterator();
		while (it.hasNext()) {
			SAMRecord r = it.next();
			if (r.isSecondaryOrSupplementary()) continue;
			for (ChimericAlignment ca : ChimericAlignment.getChimericAlignments(r)) {
				SimpleBEDFeature f = getSpanningDeletion(r, ca, minLengthPortion, minMapq);
				if (f != null) {
					writer.write(f.getContig());
					writer.write('\t');
					writer.write(Integer.toString(f.getStart() + 1));
					writer.write('\t');
					writer.write(Integer.toString(f.getEnd() + 1));
		    		writer.write('\n');
				}
			}
		}
		writer.close();
    }
    public static SimpleBEDFeature getSpanningDeletion(SAMRecord r, ChimericAlignment ca, double minLengthPortion, int minMapq) {
    	if (r.getMappingQuality() < minMapq) return null;
    	if (ca.mapq < minMapq) return null;
    	if (!ca.rname.equals(r.getReferenceName())) return null;
    	if (ca.isNegativeStrand != r.getReadNegativeStrandFlag()) return null;
    	if (r.getCigar().getReferenceLength() < r.getReadLength() * minLengthPortion) return null;
    	if (ca.cigar.getReferenceLength() < r.getReadLength() * minLengthPortion) return null;
    	int caStartClipLength = CigarUtil.getStartClipLength(ca.cigar.getCigarElements());
    	int rStartClipLength = CigarUtil.getStartClipLength(r.getCigar().getCigarElements());
    	if (r.getAlignmentEnd() < ca.pos && caStartClipLength > rStartClipLength) {
    		return new SimpleBEDFeature(r.getAlignmentEnd(), ca.pos, r.getReferenceName());
    	}
    	int endpos = ca.pos + ca.cigar.getReferenceLength() - 1;
    	if (endpos < r.getAlignmentStart() && caStartClipLength < rStartClipLength) {
    		return new SimpleBEDFeature(endpos, r.getAlignmentStart(), r.getReferenceName());
    	}
    	return null;
    }
	public static void main(String[] argv) {
	    System.exit(new SplitReadToBed().instanceMain(argv));
	}
}
