package au.edu.wehi.idsv.repeatmasker;

import au.edu.wehi.idsv.sam.CigarUtil;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.readers.LineIterator;
import joptsimple.internal.Strings;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class RepeatMaskerCodec extends AsciiFeatureCodec<RepeatMaskerFeature> {
    public static final int START_COLUMN_MATCH = 27;
    public static final int START_COLUMN_START_POSITION = 15;
    private boolean isSummaryFile = false;
    /**
     * Caching string to reduce memory usage of repeated strings
     */
    private Map<String, String> stringCache = new HashMap<>();

    public RepeatMaskerCodec() {
        super(RepeatMaskerFeature.class);
    }

    public static int parseIntIgnoringParentheses(String s) {
        return Integer.parseInt(s.substring(s.startsWith("(") ? 1 : 0, s.length() - (s.endsWith(")") ? 1 : 0)));
    }

    protected void skipEmptyLines(LineIterator lineIterator) {
        while (lineIterator.hasNext() && lineIterator.peek().equals("")) {
            lineIterator.next();
        }
    }

    @Override
    public RepeatMaskerFeature decode(LineIterator lineIterator) {
        skipEmptyLines(lineIterator);
        if (!lineIterator.hasNext()) return null;
        RepeatMaskerFeature f = new RepeatMaskerFeature();
        String line = lineIterator.next();
        skipEmptyLines(lineIterator);
        while (isHeader(line)) {
            skipEmptyLines(lineIterator);
            if (!lineIterator.hasNext()) return null;
            line = lineIterator.next();
        }
        if (!decodeSummary(line, f)) {
            skipEmptyLines(lineIterator);
            while (lineIterator.hasNext() && (lineIterator.peek().startsWith("C") || lineIterator.peek().startsWith(" "))) {
                decodeAlignment(f, lineIterator);
                skipEmptyLines(lineIterator);
            }
            while (lineIterator.hasNext() && !startsWithDigit(lineIterator.peek())) {
                lineIterator.next();
                skipEmptyLines(lineIterator);
            }
            // too brittle
//            while (lineIterator.hasNext() && (
//                    lineIterator.peek().startsWith("Matrix") ||
//                            lineIterator.peek().startsWith("Kimura") ||
//                            lineIterator.peek().startsWith("Transitions") ||
//                            lineIterator.peek().startsWith("Gap_init") ||
//                            // Footer
//                            lineIterator.peek().startsWith("##") ||
//                            lineIterator.peek().startsWith("RepeatMasker") ||
//                            lineIterator.peek().startsWith("run"))) {
//                lineIterator.next();
//                skipEmptyLines(lineIterator);
//            }
        }
        return f;
    }

    private static boolean startsWithDigit(String s) {
        if (Strings.isNullOrEmpty(s)) return false;
        return s.charAt(0) >= '0' && s.charAt(0) <= '9';
    }

    @Override
    public RepeatMaskerFeature decode(String s) {
        return decodeSummary(s);
    }

    private boolean decodeSummary(String s, RepeatMaskerFeature f) {
        String[] fields = s.trim().split(" +");
        int i = 0;
        boolean isSummary;
        f.setRepeatAlignmentSummaryInformation(new RepeatAlignmentSummaryInformation());
        f.setScore(Integer.parseInt(fields[i++]));
        f.getRepeatAlignmentSummaryInformation().setPercentageSubstituted(Float.parseFloat(fields[i++]));
        f.getRepeatAlignmentSummaryInformation().setPercentageDeleted(Float.parseFloat(fields[i++]));
        f.getRepeatAlignmentSummaryInformation().setPercentageInserted(Float.parseFloat(fields[i++]));
        f.setContig(fields[i++]);
        f.setStart(parseIntIgnoringParentheses(fields[i++]));
        f.setEnd(parseIntIgnoringParentheses(fields[i++]));
        f.getRepeatAlignmentSummaryInformation().setBasesInQueryPastMatch(parseIntIgnoringParentheses(fields[i++]));
        String strandOrRepeat = fields[i++];
        if (strandOrRepeat.equals("C")) {
            f.setStrand(Strand.NEGATIVE);
            strandOrRepeat = fields[i++];
        } else if (strandOrRepeat.equals("+")) {
            f.setStrand(Strand.POSITIVE);
            strandOrRepeat = fields[i++];
        } else {
            // align file doesn't inlcude "+", only "C"
            f.setStrand(Strand.POSITIVE);
        }
        if (strandOrRepeat.contains("#")) {
            // align format
            String[] repeatSplit = strandOrRepeat.split("#");
            f.setRepeatType(repeatSplit[0]);
            f.setRepeatClass(repeatSplit[1]);
            isSummary = false;
        } else {
            isSummary = true;
            f.setRepeatType(strandOrRepeat);
            f.setRepeatClass(fields[i++]);
        }
        if (f.getStrand() == Strand.POSITIVE) {
            f.getRepeatAlignmentSummaryInformation().setMatchStart(parseIntIgnoringParentheses(fields[i++]));
            f.getRepeatAlignmentSummaryInformation().setMatchEnd(parseIntIgnoringParentheses(fields[i++]));
            f.getRepeatAlignmentSummaryInformation().setBasesInRepeatPastMatch(parseIntIgnoringParentheses(fields[i++]));
        } else {
            f.getRepeatAlignmentSummaryInformation().setBasesInRepeatPastMatch(parseIntIgnoringParentheses(fields[i++]));
            f.getRepeatAlignmentSummaryInformation().setMatchEnd(parseIntIgnoringParentheses(fields[i++]));
            f.getRepeatAlignmentSummaryInformation().setMatchStart(parseIntIgnoringParentheses(fields[i++]));
        }
        f.setUniqueID(fields[i++]);
        return isSummary;
    }

    public RepeatMaskerFeature decodeSummary(String s) {
        RepeatMaskerFeature f = new RepeatMaskerFeature();
        decodeSummary(s, f);
        return f;
    }

    private void decodeAlignment(RepeatMaskerFeature f, LineIterator lineIterator) {
        String queryLine = lineIterator.next();
        String matchLine = lineIterator.next();
        String repeatLine = lineIterator.next();
        List<CigarElement> ce = new ArrayList<>();
        int queryStartOffset = Integer.parseInt(queryLine.substring(START_COLUMN_START_POSITION, START_COLUMN_MATCH - 1).trim());
        int queryEndOffset = Integer.parseInt(queryLine.substring(queryLine.lastIndexOf(' ') + 1));
        int repeatStartOffset = Integer.parseInt(repeatLine.substring(START_COLUMN_START_POSITION, START_COLUMN_MATCH - 1).trim());
        int repeatEndOffset = Integer.parseInt(repeatLine.substring(repeatLine.lastIndexOf(' ') + 1));
        int sequenceLength = queryLine.lastIndexOf(' ') - START_COLUMN_MATCH;
        int nestedBases = 0;
        for (int i = 0; i < sequenceLength; i++) {
            if (queryLine.charAt(START_COLUMN_MATCH + i) == 'X') {
                // This is a nested alignment
                // TODO: this code doesn't handle multiple nested Xs on the sample alignment line since we don't know
                //  how long each of the Xs are.
                int nestingLength = Math.abs(queryEndOffset - queryStartOffset) + 1 - sequenceLength + 1;
                nestedBases += nestingLength;
                ce.add(new CigarElement(nestingLength, CigarOperator.INSERTION));
            } else if (queryLine.charAt(START_COLUMN_MATCH + i) == '-') {
                ce.add(new CigarElement(1, CigarOperator.DELETION));
            } else if (repeatLine.charAt(START_COLUMN_MATCH + i) == '-') {
                ce.add(new CigarElement(1, CigarOperator.INSERTION));
            } else {
                switch (matchLine.charAt(START_COLUMN_MATCH + i)) {
                    case ' ':
                        ce.add(new CigarElement(1, CigarOperator.EQ));
                        break;
                    case 'v':
                    case 'i':
                        ce.add(new CigarElement(1, CigarOperator.X));
                        break;
                    case '?':
                        ce.add(new CigarElement(1, CigarOperator.MATCH_OR_MISMATCH));
                        break;
                    default:
                        throw new RuntimeException("Unable to parse line" + matchLine);
                }
            }
        }
        RepeatAlignmentInformation aln = f.getRepeatAlignmentInformation(false);
        if (aln == null) {
            ce.add(0, new CigarElement(queryStartOffset - 1, CigarOperator.SOFT_CLIP));
            aln = new RepeatAlignmentInformation();
            aln.setCigar(new Cigar(CigarUtil.clean(ce)));
            aln.setRepeatStart(Math.min(repeatStartOffset, repeatEndOffset));
            aln.setNestedBases(nestedBases);
            f.setAlignment(aln);
        } else {
            ce.addAll(0, aln.getCigar().getCigarElements());
            aln.setRepeatStart(Math.min(aln.getRepeatStart(), Math.min(repeatStartOffset, repeatEndOffset)));
            aln.setCigar(new Cigar(CigarUtil.clean(ce)));
            aln.setNestedBases(aln.getNestedBases() + nestedBases);
        }
    }

    @Override
    public Object readActualHeader(LineIterator lineIterator) {
        while (lineIterator.hasNext() && isHeader(lineIterator.peek())) {
            lineIterator.next();
        }
        skipEmptyLines(lineIterator);
        return null;
    }
    public static boolean isHeader(String s) {
        if (s.startsWith("There were no repetitive sequences detected")) return true;
        if (!couldBeHeader(s)) return false;
        s = s.replaceAll("\\s","");
        // Header name for first field changes depending on whether RepeatMasker is using
        // rmblast or hmmer
        boolean firstHeaderLine = s.contains("percpercpercquerypositioninquerymatchingrepeatpositioninrepeat");
        boolean secondHeaderLine = s.contains("div.del.ins.sequencebeginend(left)repeatclass/familybeginend(left)ID");
        return firstHeaderLine || secondHeaderLine;
    }

    /**
     * short-cut for quick return for non-header lines with minimal parsing
     */
    private static boolean couldBeHeader(String s) {
        for (int i = 0; i < s.length(); i++) {
            switch (s.charAt(i)) {
                case ' ':
                    continue;
                // SW
                case 'S':
                case 's':
                // bit score
                case 'b':
                case 'B':
                    return true;
                default:
                    return false;
            }
        }
        return false;
    }
    @Override
    public boolean canDecode(String s) {
        return true;
    }
}
