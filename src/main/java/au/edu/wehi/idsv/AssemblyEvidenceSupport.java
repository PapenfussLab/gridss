package au.edu.wehi.idsv;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Ordering;
import com.google.common.collect.Range;

public class AssemblyEvidenceSupport {
    public enum SupportType {
        Read(0),
        ReadPair(1);

        private final int value;

        SupportType(int value) {
            this.value = value;
        }

        public int getValue() {
            return value;
        }

        public static SupportType value(int i) {
            for (SupportType st : SupportType.values()) {
                if (st.getValue() == i) {
                    return st;
                }
            }
            throw new IllegalArgumentException("Invalid value");
        }
    }
    private final SupportType supportType;
    private final Range<Integer> assemblyContigOffset;
    private final String evidenceID;
    private final String fragmentID;
    private final int category;
    private final float qual;

    public SupportType getSupportType() {
        return supportType;
    }

    public Range<Integer> getAssemblyContigOffset() {
        return assemblyContigOffset;
    }

    public String getEvidenceID() {
        return evidenceID;
    }

    public String getFragmentID() {
        return fragmentID;
    }

    public int getCategory() {
        return category;
    }

    public AssemblyEvidenceSupport adjustForAssemblyTruncation(int startBasesTruncated) {
        return new AssemblyEvidenceSupport(
                supportType,
                Range.closed(assemblyContigOffset.lowerEndpoint() - startBasesTruncated, assemblyContigOffset.upperEndpoint() - startBasesTruncated),
                evidenceID,
                fragmentID,
                category,
                qual);
    }

    public float getQual() {
        return qual;
    }
    public AssemblyEvidenceSupport(
            SupportType supportType,
            Range<Integer> assemblyContigOffset,
            String evidenceID,
            String fragmentID,
            int category,
            float qual) {
        this.supportType = supportType;
        this.assemblyContigOffset = assemblyContigOffset;
        this.evidenceID = evidenceID;
        this.fragmentID = fragmentID;
        this.category = category;
        this.qual = qual;
    }
    public AssemblyEvidenceSupport(DirectedEvidence e, Range<Integer> supportInterval) {
        this(e instanceof NonReferenceReadPair ? SupportType.ReadPair : SupportType.Read,
                supportInterval,
                e.getEvidenceID(),
                e.getOriginatingFragmentID(((SAMEvidenceSource)(e.getEvidenceSource())).getSourceCategory()).iterator().next(),
                ((SAMEvidenceSource)(e.getEvidenceSource())).getSourceCategory(),
                e.getBreakendQual());
    }
    public static Ordering<AssemblyEvidenceSupport> ByFragmentID = new Ordering<AssemblyEvidenceSupport>() {
        public int compare(AssemblyEvidenceSupport arg1, AssemblyEvidenceSupport arg2) {
            return ComparisonChain.start()
                    .compare(arg1.fragmentID, arg2.fragmentID)
                    .result();
        }
    };
}
