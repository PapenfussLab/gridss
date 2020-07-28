package au.edu.wehi.idsv.kraken;

import au.edu.wehi.idsv.VcfBreakendSummary;
import au.edu.wehi.idsv.vcf.VcfInfoAttributes;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import java.util.Iterator;

/**
 * Annotates breakend/breakpoint sequences with NCBI taxonomic identifiers.
 * krakenInput and krakenOutput must be pipes connecting to Kraken2 input and output files specified on the command line.
 */
public class KrakenAnnotateVcf extends KrakenAnnotate<VariantContext, VariantContext> {
    private final SAMSequenceDictionary dict;

    /**
     *
     * @param krakenOutput kraken2 output parser.
     *                     The underlying stream is not closed at end of iterator.
     * @param krakenInput Fastq writer connected to kraken2 input.
     *                    Flushed if implements Flushable.
     * @param input VCF records to annotate
     */
    public KrakenAnnotateVcf(
            KrakenParser krakenOutput,
            FastqWriter krakenInput,
            Iterator<VariantContext> input,
            SAMSequenceDictionary dict,
            int minSequenceLength) {
        super(krakenOutput, krakenInput, input, minSequenceLength);
        this.dict = dict;
    }

    @Override
    protected VariantContext transform(VariantContext record, KrakenClassification kc) {
        if (kc != null && kc.isClassified) {
            record = new VariantContextBuilder(record)
                    .attribute(VcfInfoAttributes.INSERTED_SEQUENCE_NCBI_TAXONOMY_ID.attribute(), kc.taxonomyId)
                    .make();
        }
        return record;
    }

    @Override
    protected FastqRecord getKrakenFastqRecord(VariantContext vc) {
        VcfBreakendSummary be = new VcfBreakendSummary(dict, vc);
        if (be != null && be.breakpointSequence != null) {
            // TODO: does Kraken die on bad qual scores?
            return new FastqRecord(vc.getID(), be.breakpointSequence.getBytes(), null, null);
        }
        return null;
    }

    @Override
    protected boolean shouldReturnUnprocessedRecords() {
        return true;
    }
}
