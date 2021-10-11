package gridss;

import au.edu.wehi.idsv.*;
import gridss.cmdline.FullEvidenceCommandLineProgram;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.Argument;
import picard.cmdline.StandardOptionDefinitions;

import java.io.*;
import java.util.concurrent.ExecutorService;

public class DumpReadSupport extends FullEvidenceCommandLineProgram {
    private static final Log log = Log.getInstance(DumpReadSupport.class);
    public static void main(String[] argv) {
        System.exit(new DumpReadSupport().instanceMain(argv));
    }
    @Argument(shortName= StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Outputs all GRIDSS breakpoint evidence to BEDPE")
    public File OUTPUT_BEDPE;
    @Argument(doc="Outputs all GRIDSS breakend evidence to BED6")
    public File OUTPUT_BED;
    protected String[] customCommandLineValidation() {
        if (OUTPUT_BEDPE == null && OUTPUT_BED == null) {
            return new String[] {"Must specify either or both of OUTPUT_BEDPE, OUTPUT_BED"};
        }
        return super.customCommandLineValidation();
    }
    @Override
    public int doWork(ExecutorService threadpool) {
        export();
        return 0;
    }
    public void export() {
        BufferedWriter bedpe = null;
        BufferedWriter bed = null;
        try {
            if (OUTPUT_BEDPE != null) {
                bedpe = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(OUTPUT_BEDPE)));
            }
            if (OUTPUT_BED != null) {
                bed = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(OUTPUT_BEDPE)));
            }
            export(bedpe, bed);
            bedpe.close();
            bed.close();
        } catch (IOException e) {
            log.error(e);
            throw new RuntimeException(e);
        }
    }

    public void export(BufferedWriter bedpe, BufferedWriter bed) throws IOException {
        String commonFields = "\ttype\tcigar\teid";
        if (bedpe != null) {
            writeBreakpointHeaderLine(bedpe);
        }
        if (bed != null) {
            writeBreakendHeaderLine(bed);
        }
        try (CloseableIterator<DirectedEvidence> it = evidence()) {
            while (it.hasNext()) {
                DirectedEvidence e = it.next();
                if (bedpe != null && e instanceof DirectedBreakpoint) {
                    writeBreakpoint(bedpe, (DirectedBreakpoint)e);

                } if (bed != null) {
                    writeBreakend(bed, e);
                }
            }
        }
    }
    public void writeBreakpointHeaderLine(BufferedWriter writer) throws IOException {
        writer.append("#chrom1\tstart1\tend1\tchrom2\tstart2\tname\tscore\tstrand1\tstrand2");
        writeCommonHeaderFields(writer);
        writer.append("\tins\tinslen");
        writer.append("\tnompos1\tnompos2");
        writer.append("\trcigar\trlength");
        writer.newLine();
    }
    public void writeBreakpoint(BufferedWriter writer, DirectedBreakpoint e) throws IOException {
        BreakpointSummary bs = e.getBreakendSummary();
        writer.append(String.format("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%.1f\t%s\t%s",
                getReference().getSequenceDictionary().getSequence(bs.referenceIndex).getSequenceName(),
                bs.start,
                bs.end,
                getReference().getSequenceDictionary().getSequence(bs.referenceIndex2).getSequenceName(),
                bs.start2,
                bs.end2,
                e.getUnderlyingSAMRecord().getReadName(),
                e.getBreakpointQual(),
                bs.direction.toBedChar(),
                bs.direction2.toBedChar()));
        writeCommon(writer, e);
        writer.append(String.format("\t%s\t%d",
                e.getUntemplatedSequence(),
                e.getUntemplatedSequence().length()));
        writer.append(String.format("\t%d\t%d",
                bs.nominal,
                bs.nominal2));
        Cigar rcigar = null;
        if (e instanceof DiscordantReadPair) {
            rcigar = ((DiscordantReadPair) e).getNonReferenceRead().getCigar();
        } else if (e instanceof SplitReadEvidence) {
            rcigar = ((SplitReadEvidence) e).getRemoteChimericAlignment().cigar;
        } else if (e instanceof IndelEvidence) {
            rcigar = ((IndelEvidence) e).getSAMRecord().getCigar();
        }
        writer.append(String.format("\t%s\t%d",
            rcigar,
            e.getBreakendSequence().length));
        writer.newLine();
    }
    public void writeBreakendHeaderLine(BufferedWriter writer) throws IOException {
        writer.append("#chrom\tstart\tend\tname\tscore\tstrand");
        writeCommonHeaderFields(writer);
        writer.append("\tins\tinslen");
        writer.append("\tnompos");
        writer.newLine();
    }
    public void writeBreakend(BufferedWriter writer, DirectedEvidence e) throws IOException {
        BreakendSummary bs = e.getBreakendSummary();
        writer.append(String.format("%s\t%d\t%d\t%d\t%.2f\t%s",
                getReference().getSequenceDictionary().getSequence(bs.referenceIndex).getSequenceName(),
                bs.start,
                bs.end,
                e.getUnderlyingSAMRecord().getReadName(),
                e.getBreakendQual(),
                bs.direction.toBedChar()));
        writeCommon(writer, e);
        writer.append(String.format("\t%s\t%d",
                new String(e.getBreakendSequence()),
                e.getBreakendSequence().length));
        writer.append(String.format("\t%d",
                bs.nominal));
        writer.newLine();
    }
    public void writeCommonHeaderFields(BufferedWriter writer) throws IOException {
        writer.append("\ttype\teid\tcigar\tlength\tnm");
    }
    public void writeCommon(BufferedWriter writer, DirectedEvidence e) throws IOException {
        writer.append(String.format("\t%s\t%s\t%s\t%d\t%d\t",
                evidenceType(e),
                e.getEvidenceID(),
                e.getUnderlyingSAMRecord().getCigarString(),
                e.getAnchorSequence().length,
                e.getUnderlyingSAMRecord().getIntegerAttribute(SAMTag.NM.name())));

    }
    public static String evidenceType(DirectedEvidence e) {
        if (AssemblyAttributes.isAssembly(e)) {
            return "as";
        }
        if (e instanceof SplitReadEvidence) {
            return "sr";
        }
        if (e instanceof IndelEvidence) {
            return "id";
        }
        if (e instanceof SingleReadEvidence) {
            return "sc";
        }
        if (e instanceof DiscordantReadPair) {
            return "dp";
        }
        if (e instanceof NonReferenceReadPair) {
            return "um";
        }
        throw new RuntimeException("Unknown evidence type" + e.getClass().toString());
    }

    public CloseableIterator<DirectedEvidence> evidence() {
        return aes().iterator();
    }
    private AggregateEvidenceSource aes() {
        AggregateEvidenceSource aes = new AggregateEvidenceSource(getContext(), getSamEvidenceSources(), getAssemblySource(), SAMEvidenceSource.EvidenceSortOrder.EvidenceStartPosition);
        return aes;
    }
}
