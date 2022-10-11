package au.edu.wehi.idsv.sim;

import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.vcf.SvType;
import au.edu.wehi.idsv.vcf.VcfFilter;
import au.edu.wehi.idsv.vcf.VcfStructuralVariantHeaderLines;
import au.edu.wehi.idsv.vcf.VcfSvConstants;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import static au.edu.wehi.idsv.vcf.SvType.INS;

/**
 * Simple genomic rearrangement event
 */
public class SimpleEvent {
    /**
     *
     * @param type Type of rearrangement event
     * @param referenceIndex chromosome
     * @param size size of event
     * @param start 1-based start position of event
     */
    public SimpleEvent(SvType type, int referenceIndex, int size, int start, String insertedSequence, String insertedSequenceSource) {
        if (type == INS) {
            assert(insertedSequence.length() == size);
        } else {
            assert(insertedSequence.length() == 0);
        }
        this.type = type;
        this.referenceIndex = referenceIndex;
        this.start = start;
        this.size = size;
        this.insertedSequence = insertedSequence;
        this.insertedSequenceSource = insertedSequenceSource;
    }
    public final int referenceIndex;
    public final SvType type;
    /**
     * Genomic position immediately prior to event
     */
    public final int start;
    /**
     * Event size in base pairs
     */
    public final int size;
    private final String insertedSequence;
    private final String insertedSequenceSource;

    /**
     * Gets sequence in 1-based genomic coordinates
     * @param ref
     * @param start
     * @param length
     * @return
     */
    private String getSeq(ReferenceLookup ref, int start, int length) {
        //String chr = ref.getSequenceDictionary().getSequence(referenceIndex).getSequenceName();
        //String seq = new String(ref.getSubsequenceAt(chr, start, start + length - 1).getBases());
        String seq = getSubsequenceAt(ref, referenceIndex, start, length);
        assert(seq.length() == length);
        return seq;
    }
    private static String getSubsequenceAt(ReferenceLookup ref, int referenceIndex, int start, int length) {
        String chr = ref.getSequenceDictionary().getSequence(referenceIndex).getSequenceName();
        byte[] seq = ref.getSequence(chr).getBases();
        return new String(Arrays.copyOfRange(seq, start - 1, start - 1 + length), StandardCharsets.US_ASCII);
    }
    public String getVariantSeq(ReferenceLookup ref, int basesBefore, int basesAfter) {
        // nominal position isn't included
        String refSeq = getSeq(ref, start + 1, size);
        String altSeq;
        switch (type) {
            case INS:
                altSeq = insertedSequence;
                break;
            case DEL:
                altSeq = "";
                break;
            case INV:
                altSeq = SequenceUtil.reverseComplement(refSeq);
                break;
            case DUP:
                altSeq = refSeq + refSeq;
                break;
            default:
                throw new RuntimeException("NYI");
        }
        String before = getSeq(ref,start - basesBefore + 1, basesBefore);
        String after = getSeq(ref, start + 1 + getGenomicWidth(), basesAfter);
        return before + altSeq + after;
    }
    /**
     * Width of event in reference genome coordinate space
     * @return
     */
    public  int getGenomicWidth() {
        return getGenomicWidth(type, size);
    }
    /**
     * Width of event in reference genome coordinate space
     * @return
     */
    public static int getGenomicWidth(SvType type, int size) {
        switch (type) {
            case INS:
                return 0;
            case DEL:
            case INV:
            case DUP:
                return size;
            case BND:
            case CNV:
            default:
                throw new RuntimeException("Not implemented by this simulator");
        }
    }
    /**
     * Adds the VCF headers used by asVariantContextBuilder()
     * @param header
     */
    public static void addVcfHeaders(VCFHeader header) {
        header.addMetaDataLine(VcfStructuralVariantHeaderLines.END);
        header.addMetaDataLine(VcfStructuralVariantHeaderLines.SV_LENGTH);
        header.addMetaDataLine(VcfStructuralVariantHeaderLines.SV_TYPE);
        header.addMetaDataLine(VcfStructuralVariantHeaderLines.EVENT_ID);
        header.addMetaDataLine(VcfStructuralVariantHeaderLines.EVENT_TYPE);
        header.addMetaDataLine(VcfFilter.REFERENCE_ALLELE.header());
        header.addMetaDataLine(new VCFInfoHeaderLine("INS_SRC", 1, VCFHeaderLineType.String, "Location of inserted sequence"));
    }
    public String getID(ReferenceLookup ref) {
        return String.format("%s.%d.%s%d", ref.getSequenceDictionary().getSequence(referenceIndex).getSequenceName(), start, type, size);
    }
    /**
     * Convert to VCF notation
     * @return
     */
    public VariantContextBuilder asVariantContextBuilder(ReferenceLookup ref, boolean useSymbolicAllele) {
        String refActual = getSeq(ref, start, 1 + getGenomicWidth());
        String refSymbolic = refActual.substring(0, 1);
        String altSymbolic = "<" + type + ">";
        String altActual = getVariantSeq(ref, 1, 0);
        VariantContextBuilder builder = new VariantContextBuilder();
        int end = start + (type == SvType.INS ? 0 : size);
        String chr = ref.getSequenceDictionary().getSequence(referenceIndex).getSequenceName();
        builder.id(getID(ref))
                .chr(chr)
                .start(start)
                .stop(end)
                .attribute(VcfSvConstants.SV_TYPE_KEY, type)
                .attribute(VcfSvConstants.SV_LENGTH_KEY, size)
                .attribute(VcfSvConstants.EVENTTYPE_KEY, type.toString());
        if (type != INS) {
            builder.attribute(VcfSvConstants.END_KEY, end);
        }
        if (refActual.equals(altActual)) {
            builder.filter(VcfFilter.REFERENCE_ALLELE.filter());
            // htsjdk doesn't support identical alleles so we need to force it to be symbolic
            useSymbolicAllele = true;
        }
        if (useSymbolicAllele) {
            builder.alleles(refSymbolic, altSymbolic);
        } else {
            builder.alleles(refActual, altActual);
        }
        if (insertedSequenceSource != null) {
            builder.attribute("INS_SRC", insertedSequenceSource);
        }
        return builder;
    }
}
