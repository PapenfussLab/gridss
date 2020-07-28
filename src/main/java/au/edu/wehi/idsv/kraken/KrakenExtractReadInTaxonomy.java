package au.edu.wehi.idsv.kraken;

import au.edu.wehi.idsv.VcfBreakendSummary;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.vcf.VcfInfoAttributes;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.RuntimeEOFException;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import java.io.Flushable;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.concurrent.LinkedBlockingDeque;

/**
 * Annotates breakend/breakpoint sequences with NCBI taxonomic identifiers.
 * krakenInput and krakenOutput must be pipes connecting to Kraken2 input and output files specified on the command line.
 */
public class KrakenExtractReadInTaxonomy extends KrakenAnnotate<SAMRecord, String> {
    private final boolean includeSoftClips;
    private final boolean[] taxonomyLookup;

    /**
     * @param krakenOutput      kraken2 output parser.
     *                          The underlying stream is not closed at end of iterator.
     * @param krakenInput       Fastq writer connected to kraken2 input.
     *                          Flushed if implements Flushable.
     * @param input             VCF records to annotate
     * @param minSequenceLength
     */
    public KrakenExtractReadInTaxonomy(KrakenParser krakenOutput, FastqWriter krakenInput, Iterator<SAMRecord> input, int minSequenceLength, boolean includeSoftClips, boolean[] taxonomyLookup) {
        super(krakenOutput, krakenInput, input, minSequenceLength);
        this.includeSoftClips = includeSoftClips;
        this.taxonomyLookup = taxonomyLookup;
    }

    @Override
    protected String transform(SAMRecord record, KrakenClassification kc) {
        if (kc != null && kc.taxonomyId != 0 && taxonomyLookup[kc.taxonomyId]) {
            return record.getReadName();
        }
        return null;
    }

    @Override
    protected FastqRecord getKrakenFastqRecord(SAMRecord record) {
        byte[] bases = null;
        byte[] quals = null;
        if (record.getReadUnmappedFlag()) {
            bases = record.getReadBases();
            quals = record.getBaseQualities();
        } else if (includeSoftClips) {
            int startClip = SAMRecordUtil.getStartSoftClipLength(record);
            int endClip = SAMRecordUtil.getStartSoftClipLength(record);
            if (startClip > endClip) {
                // start clip
                bases = Arrays.copyOf(record.getReadBases(), startClip);
                quals = Arrays.copyOf(record.getBaseQualities(), startClip);
            } else if (endClip > 0 & endClip >= startClip) {
                // end clip
                bases = Arrays.copyOfRange(record.getReadBases(), record.getReadLength() - endClip, record.getReadLength());
                quals = Arrays.copyOfRange(record.getBaseQualities(), record.getReadLength() - endClip, record.getReadLength());
            }
        }
        if (bases != null) {
            return new FastqRecord(
                    record.getReadName(),
                    new String(bases, StandardCharsets.UTF_8),
                    null,
                    SAMUtils.phredToFastq(quals));
        }
        return null;
    }

    @Override
    protected boolean shouldReturnUnprocessedRecords() {
        return false;
    }
}
