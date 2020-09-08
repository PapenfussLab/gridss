package au.edu.wehi.idsv.repeatmasker;

import au.edu.wehi.idsv.VcfBreakendSummary;
import au.edu.wehi.idsv.sam.CigarUtil;
import au.edu.wehi.idsv.vcf.VcfInfoAttributes;
import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeaderLineType;
import org.apache.commons.compress.utils.Lists;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Annotates variants sequences with NCBI taxonomic identifiers.
 * Kraken classifications must be in the same order as the VCF records.
 */
public class AnnotateRepeatMasker /*implements Iterator<VariantContext>*/ { // can't make this an iterator since RM doesn't respect record ordering
    private static final Log log = Log.getInstance(AnnotateRepeatMasker.class);
    public static final List<VcfInfoAttributes> REPEAT_MASKER_ATTRIBUTES = Arrays.asList(
            VcfInfoAttributes.INSERTED_SEQUENCE_REPEATMASKER_ORIENTATION,
            VcfInfoAttributes.INSERTED_SEQUENCE_REPEATMASKER_REPEAT_CLASS,
            VcfInfoAttributes.INSERTED_SEQUENCE_REPEATMASKER_REPEAT_TYPE,
            VcfInfoAttributes.INSERTED_SEQUENCE_REPEATMASKER_OVERLAP,
            VcfInfoAttributes.INSERTED_SEQUENCE_REPEATMASKER_SA_TAG);
    private static final List<String> TAGS = REPEAT_MASKER_ATTRIBUTES.stream().map(a -> a.attribute()).collect(Collectors.toList());
    private final PeekingIterator<RepeatMaskerFeature> rmit;
    private final Iterator<VariantContext> it;
    private final SAMSequenceDictionary dict;
    private final List<String> tags;


    public AnnotateRepeatMasker(Iterator<RepeatMaskerFeature> rmit, Iterator<VariantContext> it, List<String> tags, SAMSequenceDictionary dict) {
        this.rmit = Iterators.peekingIterator(rmit);
        this.it = it;
        this.dict = dict;
        this.tags = tags;
    }

    public static VariantContext annotate(SAMSequenceDictionary dict, List<String> tags, VariantContext vc, Collection<RepeatMaskerFeature> rm) {
        if (rm == null || rm.isEmpty()) return vc;
        if (tags == null) tags = TAGS;
        if (rm.size() > 1) {
            // longest first
            rm = new ArrayList<>(rm);
            ((ArrayList<RepeatMaskerFeature>)rm).sort(RepeatMaskerFeature.ByAlignmentLength.reversed());
        }
        RepeatMaskerFeature f = rm.iterator().next();
        VariantContextBuilder builder = new VariantContextBuilder(vc);
        for (String tag : tags) {
            if (VcfInfoAttributes.INSERTED_SEQUENCE_REPEATMASKER_REPEAT_TYPE.attribute().equals(tag)) {
                builder = builder.attribute(VcfInfoAttributes.INSERTED_SEQUENCE_REPEATMASKER_REPEAT_TYPE.attribute(), f.getRepeatType());
            } else if (VcfInfoAttributes.INSERTED_SEQUENCE_REPEATMASKER_REPEAT_CLASS.attribute().equals(tag)) {
                builder = builder.attribute(VcfInfoAttributes.INSERTED_SEQUENCE_REPEATMASKER_REPEAT_CLASS.attribute(), f.getRepeatClass());
            } else if (VcfInfoAttributes.INSERTED_SEQUENCE_REPEATMASKER_ORIENTATION.attribute().equals(tag)) {
                builder = builder.attribute(VcfInfoAttributes.INSERTED_SEQUENCE_REPEATMASKER_ORIENTATION.attribute(), f.getStrand().toString());;
            } else if (VcfInfoAttributes.INSERTED_SEQUENCE_REPEATMASKER_OVERLAP.attribute().equals(tag)) {
                VcfBreakendSummary vbs = new VcfBreakendSummary(dict, vc);
                if (vbs.breakpointSequence.length() != 0) {
                    builder = builder.attribute(VcfInfoAttributes.INSERTED_SEQUENCE_REPEATMASKER_OVERLAP.attribute(), Math.abs(f.getEnd() - f.getStart() + 1) / (double)vbs.breakpointSequence.length());
                } else {
                    builder.rmAttribute(VcfInfoAttributes.INSERTED_SEQUENCE_REPEATMASKER_OVERLAP.attribute());
                }
            } else if (VcfInfoAttributes.INSERTED_SEQUENCE_REPEATMASKER_SA_TAG.attribute().equals(tag)) {
                builder = builder.attribute(VcfInfoAttributes.INSERTED_SEQUENCE_REPEATMASKER_SA_TAG.attribute(), toSA(rm));
            } else {
                throw new IllegalArgumentException(String.format("%s is not a known GRIDSS RepeatMasker INFO field.", tag));
            }
        }
        return builder.make();
    }

    private static List<String> toSA(Collection<RepeatMaskerFeature> rm) {
        List<String> sa = new ArrayList<>(rm.size());
        for (RepeatMaskerFeature f : rm) {
            StringBuilder sb = new StringBuilder();
            sb.append(f.getRepeatType());
            sb.append('#');
            sb.append(f.getRepeatClass());
            sb.append('|');
            sb.append(f.getRepeatAlignmentInformation(true).getRepeatStart());
            sb.append('|');
            sb.append(f.getStrand().toString());
            sb.append('|');
            sb.append(f.getRepeatAlignmentInformation(true).getCigar().toString());
            sb.append('|');
            sb.append('|');
            if (f.getRepeatAlignmentInformation(false) != null) {
                sb.append(CigarUtil.editDistance(f.getRepeatAlignmentInformation(false).getCigar(), false, false));
            }
            sa.add(sb.toString());
        }
        return sa;
    }
}
