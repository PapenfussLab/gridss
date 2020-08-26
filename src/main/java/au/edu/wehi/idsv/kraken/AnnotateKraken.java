package au.edu.wehi.idsv.kraken;

import au.edu.wehi.idsv.vcf.VcfInfoAttributes;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import joptsimple.internal.Strings;

import java.util.Iterator;

/**
 * Annotates variants sequences with NCBI taxonomic identifiers.
 * Kraken classifications must be in the same order as the VCF records.
 */
public class AnnotateKraken implements CloseableIterator<VariantContext> {
    private static final Log log = Log.getInstance(AnnotateKraken.class);
    private final KrakenParser kraken;
    private final Iterator<VariantContext> it;
    private KrakenClassification nextClassification;

    public AnnotateKraken(
            KrakenParser kraken,
            Iterator<VariantContext> it) {
        this.kraken = kraken;
        this.it = it;
    }

    @Override
    public boolean hasNext() {
        return it.hasNext();
    }

    @Override
    public VariantContext next() {
        VariantContext vc = it.next();
        return annotate(vc);
    }

    private VariantContext annotate(VariantContext vc) {
        while (nextClassification == null && kraken.hasNext()) {
            nextClassification = kraken.next();
        }
        if (nextClassification != null && nextClassification.sequenceId.equals(vc.getID())) {
            if (Strings.isNullOrEmpty(nextClassification.sequenceId) || nextClassification.sequenceId.equals(".")) {
                log.error("Kraken classification record has identifier of '" + nextClassification.sequenceId + "'.");
            }
            vc = annotate(vc, nextClassification);
            nextClassification = null;
        }
        return vc;
    }

    private static VariantContext annotate(VariantContext vc, KrakenClassification classification) {
        if (classification == null || !classification.isClassified) {
            return vc;
        }
        vc = new VariantContextBuilder(vc)
                .attribute(VcfInfoAttributes.INSERTED_SEQUENCE_NCBI_TAXONOMY_ID.attribute(), classification.taxonomyId)
                .make();
        return vc;
    }

    @Override
    public void close() {
        CloserUtil.close(it);
        CloserUtil.close(kraken);
    }
}
