package au.edu.wehi.idsv.kraken;

import au.edu.wehi.idsv.ncbi.MinimalTaxonomyNode;
import au.edu.wehi.idsv.ncbi.TaxonomyHelper;
import gridss.kraken.SubsetToTaxonomy;
import htsjdk.samtools.util.Log;
import org.apache.commons.collections4.iterators.ReverseListIterator;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

/**
 * Checks whether this read could be part of the given taxonomic subset.
 *
 * This method check not only the overall classification but also checks
 * whether this read represent a split read between the taxonomic subset
 * and an outgroup.
 *
 * This can be used to identify reads spanning a viral integration site.
 */
public class KrakenClassificationChecker {
    private static final Log log = Log.getInstance(KrakenClassificationChecker.class);
    /**
     * Taxa that are definitely in our set of interest.
     * Ancestors of taxa in our set of interest are ambiguous.
     */
    private final boolean[] goodTaxId;
    /**
     * Taxa that are definitely not in our set of interest.
     * Ancestors of taxa in our set of interest are ambiguous.
     */
    private final boolean[] badTaxId;
    public KrakenClassificationChecker(List<Integer> taxonomyIdOfInterest, File nodesdmp) throws IOException {
        log.info("Loading NCBI taxonomy from ", nodesdmp);
        Map<Integer, MinimalTaxonomyNode> lookup = TaxonomyHelper.parseMinimal(nodesdmp);
        this.goodTaxId = TaxonomyHelper.createInclusionLookup(taxonomyIdOfInterest, lookup);
        this.badTaxId = setupBadTaxId(goodTaxId, lookup, taxonomyIdOfInterest);
    }

    private boolean[] setupBadTaxId(boolean[] goodTaxId, Map<Integer, MinimalTaxonomyNode> lookup, List<Integer> taxonomyIdOfInterest) {
        boolean[] badTaxId = new boolean[goodTaxId.length];
        for (int i = 0 ; i < goodTaxId.length; i++) {
            badTaxId[i] = !goodTaxId[i];
        }
        for (int ancestorTaxId : taxonomyIdOfInterest) {
            int lastTaxId = -1;
            while (ancestorTaxId != lastTaxId) {
                badTaxId[ancestorTaxId] = false;
                lastTaxId = ancestorTaxId;
                ancestorTaxId = lookup.get(ancestorTaxId).parentTaxId;

            }
        }
        return badTaxId;
    }

    public boolean isOfInterest(KrakenClassification kc) {
        if (goodTaxId[kc.taxonomyId]) return true;
        return isOfInterest(kc.kmerTaxonomyIds) || isOfInterest(kc.kmerTaxonomyIds2);
    }
    private boolean isOfInterest(List<KrakenKmerClassification> read) {
        return !read.isEmpty() && (
                isOfInterest(read.iterator()) || isOfInterest(new ReverseListIterator<>(read)));
    }
    private boolean isOfInterest(Iterator<KrakenKmerClassification> it) {
        // Traverse towards the middle of the read
        // if we find at least one good kmer and any number of ambiguous kmers
        // then we might be a split read to a taxonomic sequence of interest
        boolean foundGood = false;
        while (it.hasNext()) {
            KrakenKmerClassification kkc = it.next();
            if (kkc.taxonomyId != KrakenKmerClassification.AMBIGUOUS) {
                foundGood |= goodTaxId[kkc.taxonomyId];
                if (badTaxId[kkc.taxonomyId]) {
                    break;
                }
            }
        }
        return foundGood;
    }
}
