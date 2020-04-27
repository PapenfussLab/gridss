package au.edu.wehi.idsv.alignment;

import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import com.google.common.collect.ImmutableList;
import htsjdk.samtools.*;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SequenceUtil;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.net.URI;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Wrapper for the BwaMemAligner that returns SAMRecords
 */
public class BwaAligner implements Closeable {
    private static final Log log = Log.getInstance(BwaAligner.class);
    private final BwaMemIndex index;
    private final BwaMemAligner aligner;
    private final SAMSequenceDictionary dict;
    private final SAMFileHeader header;

    public BwaMemAligner getAligner() {
        return this.aligner;
    }

    public BwaAligner(File reference, SAMSequenceDictionary dict, int threads) {
        this.index = getBwaIndexFor(reference);
        this.dict = dict;
        this.header = getMinimalHeader(dict);
        this.aligner = new BwaMemAligner(this.index);
        this.aligner.setNThreadsOption(threads);
        this.aligner.setClip3PenaltyOption(0);
        this.aligner.setClip5PenaltyOption(0);
        try {
            ensureMatchingReferences(this.index, dict);
        } catch (IllegalArgumentException e) {
            // don't leak the index since it's huge
            close();
            throw e;
        }
    }

    public static BwaMemIndex getBwaIndexFor(File reference) {
        File image = new File(reference.getAbsolutePath() + BwaMemIndex.IMAGE_FILE_EXTENSION);
        if (!image.exists()) {
            log.warn("Unable to find " + image.toString() + ". Attempting to create.");
            if (BwaMemIndex.INDEX_FILE_EXTENSIONS.stream().allMatch(suffix -> new File(reference.getAbsolutePath() + suffix).exists())) {
                log.warn("Found bwa index files. Attempting to create image from bwa index files");
                System.err.flush();
                BwaMemIndex.createIndexImageFromIndexFiles(reference.getAbsolutePath(), image.getAbsolutePath());
            } else {
                log.warn("Could not find bwa index files. Creating bwa image from reference genome. This is a one-time operation and may take several hours.");
                System.err.flush();
                BwaMemIndex.createIndexImageFromFastaFile(reference.getAbsolutePath(), image.getAbsolutePath());
            }
            if (image.exists()) {
                log.info("Index creation successful");
            } else {
                String msg = "Index creation failed for index file " + image.toString();
                log.error(msg);
                throw new RuntimeException(msg);
            }
        }
        log.info("Loading bwa mem index image from " + image);
        System.err.flush(); // ensure our warning error message gets to the console as we're possible about to die in C code
        BwaMemIndex index = new BwaMemIndex(image.getAbsolutePath());
        return index;
    }

    private static SAMFileHeader getMinimalHeader(SAMSequenceDictionary dict) {
        SAMFileHeader header = new SAMFileHeader();
        for (SAMSequenceRecord ref : dict.getSequences()) {
            header.addSequence(ref);
        }
        return header;
    }

    public static void ensureMatchingReferences(BwaMemIndex index, SAMSequenceDictionary dict) {
        String indexNames = index.getReferenceContigNames().stream().collect(Collectors.joining("   "));
        String refNames = dict.getSequences().stream().map(x -> x.getSequenceName()).collect(Collectors.joining("   "));
        if (!indexNames.equals(refNames)) {
            throw new IllegalArgumentException("bwa index and reference genome sequences do not match");
        }
    }

    public List<SAMRecord> align(Collection<FastqRecord> input) {
        List<byte[]> inputs = new ArrayList<>(input.size());
        for (FastqRecord fq : input) {
            inputs.add(fq.getReadBases());
        }
        log.debug(String.format("Aligning %d sequences using BWA JNI", inputs.size()));
        List<List<BwaMemAlignment>> bwaResult = aligner.alignSeqs(inputs);
        if (bwaResult.size() != input.size()) {
            throw new IllegalStateException(String.format("bwa returned alignments for %d reads, when input with %d reads.", bwaResult.size(), input.size()));
        }
        List<SAMRecord> samResult = new ArrayList<>((int)(input.size() * 1.3)); // conservatively guess 30% of alignments are split read alignments
        int i = 0;
        for (FastqRecord fq : input) {
            List<BwaMemAlignment> bma = bwaResult.get(i++);
            List<SAMRecord> alignments = transform(fq, bma);
            samResult.addAll(alignments);
        }
        return samResult;
    }

    public List<SAMRecord> transform(FastqRecord fq, List<BwaMemAlignment> bma) {
        List<SAMRecord> result = new ArrayList<>(bma.size() == 0 ? 1 : bma.size());
        if (bma.size() == 0 || bma.get(0).getRefId() == -1) {
            SAMRecord r = SAMRecordUtil.createSAMRecord(header, fq, false);
            r.setReadUnmappedFlag(true);
            result.add(r);
        } else {
            for (BwaMemAlignment alignment : bma) {
                SAMRecord r = createAlignment(fq, alignment);
                result.add(r);
            }
            SAMRecordUtil.reinterpretAsSplitReadAlignment(result, true);
        }
        return  result;
    }

    private SAMRecord createAlignment(FastqRecord fq, BwaMemAlignment alignment) {
        SAMRecord r = SAMRecordUtil.createSAMRecord(header, fq, (alignment.getSamFlag() & SAMFlag.READ_REVERSE_STRAND.intValue()) != 0);

        r.setFlags(alignment.getSamFlag());
        r.setReferenceIndex(alignment.getRefId());
        r.setAlignmentStart(alignment.getRefStart() + 1);
        r.setCigarString(alignment.getCigar());
        if (r.getCigar().getReadLength() != fq.getReadLength()) {
            // TODO: do I need to add soft clip padding to the alignment cigar?
            int leftClipping = alignment.getSeqStart();
            int rightClipping = fq.getReadLength() - alignment.getSeqEnd();
            //r.setCigarString( CigarUtil.addSoftClipping(leftClipping, rightClipping));
            throw new IllegalStateException(String.format("Read length is %d, cigar is %s", fq.getReadLength(), r.getCigarString()));
        }
        r.setMappingQuality(alignment.getMapQual());
        r.setMateReferenceIndex(alignment.getMateRefId());
        r.setMateReferenceIndex(alignment.getMateRefStart() + 1);
        r.setInferredInsertSize(alignment.getTemplateLen());
        r.setAttribute(SAMTag.MD.name(), alignment.getMDTag());
        r.setAttribute(SAMTag.NM.name(), alignment.getNMismatches());
        r.setAttribute(SAMTag.AS.name(), alignment.getAlignerScore());
        r.setAttribute("XA", alignment.getXATag());
        r.setAttribute("XS", alignment.getSuboptimalScore());
        return r;
    }

    @Override
    public void close() {
        this.index.close();
    }
}
