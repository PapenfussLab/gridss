package au.edu.wehi.idsv;

import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.Log;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public abstract class SplitReadRealigner {
    private static final Log log = Log.getInstance(SplitReadRealigner.class);
    private final ReferenceLookup reference;
    private EvidenceIdentifierGenerator eidgen = new HashedEvidenceIdentifierGenerator();
    private int workerThreads = Runtime.getRuntime().availableProcessors();
    private int minSoftClipLength = 1;
    private float minSoftClipQuality = 0;
    private boolean realignExistingSplitReads = false;
    private boolean realignEntireRecord = false;
    private boolean processSecondaryAlignments = false;
    private boolean realignAnchoringBases = false;
    private boolean adjustPrimary = false;
    private boolean writeOA = true;

    public SplitReadRealigner(ReferenceLookup reference) {
        this.reference = reference;
    }

    public abstract void createSupplementaryAlignments(File input, File output, File unorderedOutput) throws IOException;

    protected boolean shouldDropInputRecord(SAMRecord r) {
        if ((isRealignEntireRecord() || isRealignExistingSplitReads()) && r.getSupplementaryAlignmentFlag()) {
            // drop existing supp alignments
            return true;
        }
        return false;
    }

    protected void writeCompletedAlignment(SAMRecord primary, List<SAMRecord> realignments, SAMFileWriter coordinateSortedWriter, SAMFileWriter unorderedWriter) {
        if (isRealignExistingSplitReads() || isRealignEntireRecord()) {
            if (primary.getSupplementaryAlignmentFlag()) {
                // If we're realigning, we need to drop all existing supplementary alignments
                return;
            }
            primary.setAttribute(SAMTag.SA.name(), null);
        }
        boolean primaryHasMoved = prepareRecordsForWriting(primary, realignments);
        if (primary.getReadUnmappedFlag() || primaryHasMoved) {
            // we'll break sort ordering if we write it back to the input file
            unorderedWriter.addAlignment(primary);
        } else {
            coordinateSortedWriter.addAlignment(primary);
        }
        for (SAMRecord sar : realignments) {
            unorderedWriter.addAlignment(sar);
        }
    }

    protected boolean prepareRecordsForWriting(SAMRecord primary, List<SAMRecord> realignments) {
        int primaryReferenceIndex = primary.getReadUnmappedFlag() ? -1 : primary.getReferenceIndex();
        int primaryAlignmentStart = primary.getReadUnmappedFlag() ? -1 : primary.getAlignmentStart();
        // Only do anything to records that we actually attempted realignment for
        if (realignments.size() > 0) {
            if (isRealignEntireRecord() &&
                    // special case exclusion of unanchored assemblies as
                    // we repurpose the CIGAR string to encode the breakend interval
                    // and realigning will break that
                    !AssemblyAttributes.isUnanchored(primary)) {
                SAMRecord newPrimaryAlignmentPosition = SplitReadHelper.replaceAlignment(primary, realignments, isWriteOATag());
                if (newPrimaryAlignmentPosition != null && !realignments.remove(newPrimaryAlignmentPosition)) {
                    throw new RuntimeException("Sanity check failure: no supplementary alignment was removed when replacing alignment");
                }
            } else if (isRealignAnchoringBases()) {
                // pull out the anchoring base alignment from the list and process it.
                for (int i = 0; i < realignments.size(); i++) {
                    SAMRecord r = realignments.get(i);
                    if (SplitReadHelper.isAnchoringBasesRecord(r)) {
                        if (!r.getSupplementaryAlignmentFlag()) {
                            SplitReadHelper.rewriteAnchor(primary, r);
                        }
                        realignments.remove(i);
                    }
                }
            }
            SplitReadHelper.convertToSplitRead(primary, realignments, getReference(), isAdjustPrimaryAlignment() || isRealignEntireRecord());
            for (int i = 0; i < realignments.size(); i++) {
                SAMRecord r = realignments.get(i);
                if (r.getReadUnmappedFlag() || SAMRecordUtil.getStartClipLength(r) == SAMRecordUtil.getReadLengthIncludingHardClipping(r)) {
                    realignments.remove(i);
                    i--;
                }
            }
            // TODO: remove alignments which are contained by another alignment
            if (primary.getReadUnmappedFlag() || SAMRecordUtil.getStartClipLength(primary) == SAMRecordUtil.getReadLengthIncludingHardClipping(primary)) {
                SplitReadHelper.replaceAlignment(primary, realignments, writeOA);
            }
        }
        return primary.getReadUnmappedFlag() || primary.getReferenceIndex() != primaryReferenceIndex || primary.getAlignmentStart() != primaryAlignmentStart;
    }

    public List<FastqRecord> extract(SAMRecord r, boolean isRecursiveRealignment) {
        if (!isRecursiveRealignment && shouldDropInputRecord(r)) {
            throw new IllegalArgumentException("Record should have been dropped.");
        }
        if (isRealignEntireRecord() && isRealignAnchoringBases()) {
            throw new IllegalArgumentException("Cannot realign anchoring bases if realigning entire read");
        }
        List<FastqRecord> list = new ArrayList<>(2);
        if (r.getReadUnmappedFlag()) return list;
        if (!SAMRecordUtil.isSoftClipLengthAtLeast(r, getMinSoftClipLength())) return list;
        if (isRealignEntireRecord() && !isRecursiveRealignment && !AssemblyAttributes.isUnanchored(r)) {
            list.add(SplitReadHelper.getFullRealignment(r, getEvidenceIdentifierGenerator()));
            return list;
        }
        if (isRecursiveRealignment && SplitReadHelper.isAnchoringBasesRecord(r)) {
            // don't chain realigning the anchoring bases alignments record - it's a once off
            return list;
        }
        // Logic for extending an existing SA alignment not yet complete. Need to:
        // - only realign bases not in any existing SA alignment
        // - update all SA record (requires queryname sorted input file)
        // Note that realignments may be split reads, but we currently only consider the primary alignment
        if (!isRealignExistingSplitReads() && !isRecursiveRealignment && r.getAttribute(SAMTag.SA.name()) != null) return list;
        if (r.getSupplementaryAlignmentFlag()) return list;
        if (r.isSecondaryAlignment() && !isProcessSecondaryAlignments()) return list;
        for (FastqRecord fqr : SplitReadHelper.getSplitReadRealignments(r, isRecursiveRealignment, getEvidenceIdentifierGenerator())) {
            if (fqr.getReadLength() < getMinSoftClipLength()) continue;
            if (averageBaseQuality(fqr) < getMinSoftClipQuality()) continue;
            list.add(fqr);
        }
        if (isRealignAnchoringBases() && !isRecursiveRealignment && !AssemblyAttributes.isUnanchored(r) && list.size() > 0) {
            // only emit the anchoring bases:
            // - on the first pass where the anchoring bases are actually the read anchoring bases
            // - when we actually have anchoring bases
            // - when we have something that we might actually realign (TODO: is this the best approach?)
            list.add(SplitReadHelper.getAnchoringBases(r, getEvidenceIdentifierGenerator()));
        }
        return list;
    }
    private static double averageBaseQuality(FastqRecord fqr) {
        long sum = 0;
        for (byte v : SAMUtils.fastqToPhred(fqr.getBaseQualityString())) {
            sum += v;
        }
        return (double)sum / fqr.getBaseQualityString().length();
    }

    public int getMinSoftClipLength() {
        return minSoftClipLength;
    }

    public void setMinSoftClipLength(int minSoftClipLength) {
        this.minSoftClipLength = minSoftClipLength;
    }

    public float getMinSoftClipQuality() {
        return minSoftClipQuality;
    }

    public void setMinSoftClipQuality(float minSoftClipQuality) {
        this.minSoftClipQuality = minSoftClipQuality;
    }

    public int getWorkerThreads() {
        return workerThreads;
    }

    public void setWorkerThreads(int workerThreads) {
        this.workerThreads = workerThreads;
    }

    public boolean isProcessSecondaryAlignments() {
        return processSecondaryAlignments;
    }

    public void setProcessSecondaryAlignments(boolean processSecondaryAlignments) {
        this.processSecondaryAlignments = processSecondaryAlignments;
    }

    public ReferenceLookup getReference() {
        return reference;
    }

    public boolean isRealignExistingSplitReads() {
        return realignExistingSplitReads;
    }

    public void setRealignExistingSplitReads(boolean realignExistingSplitReads) {
        this.realignExistingSplitReads = realignExistingSplitReads;
    }

    public boolean isRealignEntireRecord() {
        return realignEntireRecord;
    }

    public void setRealignEntireRecord(boolean realignEntireRecord) {
        this.realignEntireRecord = realignEntireRecord;
    }

    public boolean isRealignAnchoringBases() { return this.realignAnchoringBases; }

    public void setRealignAnchoringBases(boolean realignAnchoringBases) {
        this.realignAnchoringBases = realignAnchoringBases;
    }

    public void setAdjustPrimaryAlignment(boolean adjustPrimary) {
        this.adjustPrimary = adjustPrimary;
    }

    public boolean isAdjustPrimaryAlignment() {
        return adjustPrimary;
    }

    /**
     * Alignment-unique read identifier must be hashed to ensure that the read names
     * written to the fastq files do not exceed the BAM limit of 254 even when the
     * input reads names are at this limit.
     * See https://github.com/PapenfussLab/gridss/issues/82
     */
    public EvidenceIdentifierGenerator getEvidenceIdentifierGenerator() {
        return eidgen;
    }

    public void setEvidenceIdentifierGenerator(EvidenceIdentifierGenerator eidgen) {
        this.eidgen = eidgen;
    }

    public boolean isWriteOATag() {
        return writeOA;
    }

    public void setWriteOATag(boolean writeOA) {
        this.writeOA = writeOA;
    }
}
