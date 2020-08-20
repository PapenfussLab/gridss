package au.edu.wehi.idsv.ncbi;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.Collection;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * In-memory representation of the NCBI Taxonomy
 */
public abstract class TaxonomyHelper {
    public static final int MAX_NCBI_TAXID = 2758539;
    /**
     * Parses nodes.dmp from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
     *
     * @param nodesdmp nodes.dmp file
     */
    public static Map<Integer, TaxonomyNode> parseFull(File nodesdmp) throws IOException {
        return Files.lines(nodesdmp.toPath())
                .map(s -> new TaxonomyNode(s))
                .collect(Collectors.toMap(tn -> tn.taxId, Function.identity()));
    }

    public static Map<Integer, MinimalTaxonomyNode> parseMinimal(File nodesdmp) throws IOException {
        return Files.lines(nodesdmp.toPath())
                .map(s -> new MinimalTaxonomyNode(s))
                .collect(Collectors.toMap(tn -> tn.taxId, Function.identity()));
    }

    public static TaxonomyNode[] toArray(Map<Integer, TaxonomyNode> lookup) {
        TaxonomyNode[] array = new TaxonomyNode[maxTaxId(lookup) + 1];
        lookup.values().stream().forEach(tn -> array[tn.taxId] = tn);
        return array;
    }

    public static int maxTaxId(Map<Integer, ? extends MinimalTaxonomyNode> lookup) {
        return lookup.values().stream().mapToInt(tn -> tn.taxId).max().orElse(0);
    }

    /**
     * Creates a lookup table indicating whether that node is included, or is a child of
     * any of the given NCBI taxonomy IDs.
     * @param taxIds taxonomy IDs to search for
     * @return lookup table of inclusion(true) or exclusion(false) any of the given taxonomy IDs.
     */
    public static boolean[] createInclusionLookup(Collection<Integer> taxIds, Map<Integer, ? extends MinimalTaxonomyNode> lookup) {
        boolean[] result = new boolean[maxTaxId(lookup) + 1];
        for (int id : taxIds) {
            result[id] = true;
        }
        // check all children
        for (int i = 0; i < result.length; i++) {
            MinimalTaxonomyNode node = lookup.get(i);
            while (node != null && node.taxId != node.parentTaxId) {
                if (result[node.taxId]) {
                    result[i] = true;
                    break;
                }
                node = lookup.get(node.parentTaxId);
            }
        }
        return result;
    }

    public static boolean[] leafNodes(Map<Integer, MinimalTaxonomyNode> lookup) {
        boolean[] result = new boolean[maxTaxId(lookup) + 1];
        Arrays.fill(result, true);
        for (MinimalTaxonomyNode n : lookup.values()) {
            result[n.parentTaxId] = false;
        }
        return result;
    }
    public static boolean[] addAncestors(boolean[] taxa, Map<Integer, MinimalTaxonomyNode> lookup) {
        boolean[] result = Arrays.copyOf(taxa, taxa.length);
        for (int i = 0; i < result.length; i++) {
            if (taxa[i]) {
                // iterate up the tree
                MinimalTaxonomyNode node = lookup.get(i);
                while (node != null && node.taxId != node.parentTaxId) {
                    node = lookup.get(node.parentTaxId);
                    result[node.taxId] = true;
                }
            }
        }
        return result;
    }
}
