package au.edu.wehi.idsv;

import au.edu.wehi.idsv.bed.RepeatMaskerBEDCodec;
import au.edu.wehi.idsv.bed.RepeatMaskerBEDFeature;
import au.edu.wehi.idsv.sam.ChimericAlignment;
import au.edu.wehi.idsv.util.IntervalUtil;
import au.edu.wehi.idsv.vcf.VcfInfoAttributes;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.Log;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

public class VariantContextRepeatMaskerAnnotator implements Function<VariantContext, VariantContext> {
    private static final Log log = Log.getInstance(VariantContextRepeatMaskerAnnotator.class);
    private Map<String, IntervalTree<RepeatMaskerBEDFeature>> lookup;

    private class RepeatMaskerHit {
        public final RepeatMaskerBEDFeature rme;
        public final double overlap;
        public final boolean isNegative;
        public RepeatMaskerHit(RepeatMaskerBEDFeature rme, int alignmentStart, int alignmentEnd, boolean alignmentOnNegative) {
            this.rme = rme;
            this.overlap = IntervalUtil.overlapsWidthClosed(alignmentStart, alignmentEnd, rme.getStart(), rme.getEnd()) / (alignmentEnd - alignmentStart + 1.0);
            this.isNegative = alignmentOnNegative != (rme.getStrand() == Strand.NEGATIVE);
        }
    }

    public Collection<String> getRepeatMaskerContigs() {
        return lookup.keySet();
    }

    private static Map<String, IntervalTree<RepeatMaskerBEDFeature>> createLookup(File repeatMaskerBed) throws IOException {
        Map<String, IntervalTree<RepeatMaskerBEDFeature>> lookup = new HashMap<>();
        try (AbstractFeatureReader<BEDFeature, LineIterator> reader = AbstractFeatureReader.getFeatureReader(repeatMaskerBed.getPath(), new RepeatMaskerBEDCodec(), false)) {
            int lineno = 0;
            for (BEDFeature rawfeat : reader.iterator()) {
                RepeatMaskerBEDFeature feat = (RepeatMaskerBEDFeature)rawfeat;
                lineno++;
                IntervalTree<RepeatMaskerBEDFeature> tree = lookup.get(feat.getContig());
                if (tree == null) {
                    tree = new IntervalTree<>();
                    lookup.put(feat.getContig(), tree);
                }
                tree.put(rawfeat.getStart(), rawfeat.getEnd(), feat);
            }
        }
        return lookup;
    }

    public VariantContextRepeatMaskerAnnotator(File repeatMaskerBed) throws IOException {
        this.lookup = createLookup(repeatMaskerBed);
    }

    private Stream<RepeatMaskerHit> getHits(String s) {
        ChimericAlignment aln = ChimericAlignment.parseBEALNAlignment(s);
        String chr = aln.rname;
        int start = aln.pos;
        int end = aln.pos + aln.cigar.getReferenceLength() - 1;
        IntervalTree<RepeatMaskerBEDFeature> tree = lookup.get(chr);
        if (tree == null) return Stream.empty();
        return StreamSupport.stream(Spliterators.spliteratorUnknownSize(tree.overlappers(start, end), 0), false)
                .map(n -> new RepeatMaskerHit(n.getValue(), start, end, aln.isNegativeStrand));
    }

    private VariantContext annotate(VariantContext variantContext, RepeatMaskerHit hit) {
        return new VariantContextBuilder(variantContext)
                .attribute(VcfInfoAttributes.INSERTED_SEQUENCE_REPEATMASKER_OVERLAP.attribute(), hit.overlap)
                .attribute(VcfInfoAttributes.INSERTED_SEQUENCE_REPEATMASKER_REPEAT_TYPE.attribute(), hit.rme.getRepeatType())
                .attribute(VcfInfoAttributes.INSERTED_SEQUENCE_REPEATMASKER_REPEAT_CLASS.attribute(), hit.rme.getRepeatClass())
                .attribute(VcfInfoAttributes.INSERTED_SEQUENCE_REPEATMASKER_ORIENTATION.attribute(), hit.isNegative ? "-" : "+")
                .make();
    }

    @Override
    public VariantContext apply(VariantContext variantContext) {
        try {
            List<String> alignments = variantContext.getAttributeAsStringList(VcfInfoAttributes.BREAKEND_ALIGNMENTS.attribute(), null);
            if (alignments == null || alignments.size() == 0) {
                return variantContext;
            }
            Optional<RepeatMaskerHit> bestHit = alignments.stream()
                    .flatMap(a -> getHits(a))
                    .sorted(ByOverlapDesc)
                    .findFirst();
            if (bestHit.isPresent()) {
                variantContext = annotate(variantContext, bestHit.get());
            }
        } catch (IndexOutOfBoundsException ioobe) {
            log.error("Malformed BEALN field for " + variantContext.getID() + ". Ignoring.");
        } catch (NumberFormatException nfe) {
            log.error("Malformed BEALN field for " + variantContext.getID() + ". Ignoring.");
        }
        return variantContext;
    }
    public static final Comparator<RepeatMaskerHit> ByOverlapDesc = Comparator.comparingDouble(hit -> -hit.overlap);
}
