package au.edu.wehi.idsv.vcf;

import au.edu.wehi.idsv.alignment.StreamingAligner;
import au.edu.wehi.idsv.alignment.StreamingAlignerIterator;
import au.edu.wehi.idsv.util.AutoClosingIterator;
import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import joptsimple.internal.Strings;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class InsertedSequenceAnnotator implements CloseableIterator<VariantContext> {
	public static final byte DEFAULT_QUAL_SCORE = 20;
	private static final Log log = Log.getInstance(InsertedSequenceAnnotator.class);
	private static final Pattern breakendRegex = Pattern.compile("^(.(?<leftins>.*))?[\\[\\]].*[\\[\\]]((?<rightins>.*).)?$");
	private final File vcf;
	private final boolean stripExistingBEALN;
	private final boolean skipExistingBEALN;
	private CloseableIterator<VariantContext> vcfStream;
	private PeekingIterator<SAMRecord> alignerStream;
	private Thread feedingAligner;
	private VariantContext nextRecord = null;
	public InsertedSequenceAnnotator(File vcf, StreamingAligner aligner, boolean stripExistingBEALN, boolean skipExistingBEALN) {
		this.vcf = vcf;
		this.stripExistingBEALN = stripExistingBEALN;
		this.skipExistingBEALN = skipExistingBEALN;
		this.vcfStream = getVcf();
		StreamingAlignerIterator sai = new StreamingAlignerIterator(aligner);
		this.alignerStream = Iterators.peekingIterator(sai);
		this.feedingAligner = new Thread(() -> feedStreamingAligner(sai, aligner));
		this.feedingAligner.setName("feedAligner");
		this.feedingAligner.start();
	}

	private CloseableIterator<VariantContext> getVcf() {
		VCFFileReader vcfReader = new VCFFileReader(vcf, false);
		CloseableIterator<VariantContext> it = vcfReader.iterator();
		return new AutoClosingIterator<VariantContext>(it, vcfReader);
	}
	private static String getBreakendSequence(VariantContext seq) {
		if (seq.getAlternateAlleles().size() != 1) return null;
		Allele allele = seq.getAlternateAllele(0);
		String alt = allele.getDisplayString();
		if (alt.charAt(0) == '.' || alt.charAt(alt.length() - 1) == '.') {
			return alt.substring(1, alt.length() - 1);
		} else if (alt.charAt(0) == '[' || alt.charAt(0) == ']' || alt.charAt(alt.length() - 1) == '[' || alt.charAt(alt.length() - 1) == ']') {
			Matcher matcher = breakendRegex.matcher(alt);
			if (matcher.matches()) {
				String leftIns = matcher.group("leftins");
				if (!Strings.isNullOrEmpty(leftIns)) {
					return leftIns;
				}
				String rightIns = matcher.group("rightins");
				return rightIns;
			}
		}
		return null;
	}
	private VariantContext stripIfNeeded(VariantContext vc) {
		if (stripExistingBEALN && vc.hasAttribute(VcfInfoAttributes.BREAKEND_ALIGNMENTS.attribute())) {
			return new VariantContextBuilder(vc)
					.rmAttribute(VcfInfoAttributes.BREAKEND_ALIGNMENTS.attribute())
					.make();
		}
		return vc;
	}
	private boolean shouldSkipRecord(VariantContext vc) {
		return skipExistingBEALN && vc.hasAttribute(VcfInfoAttributes.BREAKEND_ALIGNMENTS.attribute());
	}
	private void feedStreamingAligner(StreamingAlignerIterator wrapper, StreamingAligner aligner) {
		try (CloseableIterator<VariantContext> it = getVcf()) {
			while (it.hasNext()) {
				VariantContext vc = stripIfNeeded(it.next());
				if (shouldSkipRecord(vc)) {
					// skip this record
				} else {
					String seqstr = getBreakendSequence(vc);
					if (!Strings.isNullOrEmpty(seqstr)) {
						byte[] seq = seqstr.getBytes(StandardCharsets.UTF_8);
						byte[] qual = new byte[seq.length];
						Arrays.fill(qual, DEFAULT_QUAL_SCORE);
						FastqRecord fq = new FastqRecord(vc.getID(), seq, null, qual);
						aligner.asyncAlign(fq);
					}
				}
			}
		} catch (IOException e) {
			log.warn(e);
		} finally {
			try {
				wrapper.flush();
				wrapper.close();
				aligner.close();
			} catch (IOException e) {
				log.warn(e);
			}
		}
		log.debug("Completed async external alignment feeder thread.");
	}

	@Override
	public boolean hasNext() {
		ensureNext();
		return nextRecord != null;
	}
	@Override
	public VariantContext next() {
		ensureNext();
		if (!hasNext()) throw new NoSuchElementException();
		VariantContext result = nextRecord;
		nextRecord = null;
		return result;
	}
	private void ensureNext() {
		if (nextRecord != null) {
			return;
		}
		if (!vcfStream.hasNext()) {
			nextRecord = null;
			if (alignerStream.hasNext()) {
				log.debug("Traversing aligner stream to enable graceful termination.");
				while(alignerStream.hasNext()) {
					// consume the aligner output so everything closes gracefully
				}
			}
			return;
		}
		nextRecord = stripIfNeeded(vcfStream.next());
		if (!shouldSkipRecord(nextRecord)) {
			annotateNextRecord();
		}
	}
	private void annotateNextRecord() {
		String readName = nextRecord.getID();
		List<SAMRecord> alignments = new ArrayList<>();
		while (alignerStream.hasNext() && alignerStream.peek().getReadName().equals(readName)) {
			SAMRecord r = alignerStream.next();
			if (!r.getReadUnmappedFlag()) {
				alignments.add(r);
			}
		}
		VariantContextBuilder builder = new VariantContextBuilder(nextRecord);
		List<String> existingAlignments = nextRecord.getAttributeAsStringList(VcfInfoAttributes.BREAKEND_ALIGNMENTS.attribute(), null);
		List<String> newAlignments = writeAlignmentAnnotation(alignments);
		List<String> mergedAlignments = new ArrayList<>();
		if (existingAlignments != null && existingAlignments.size() > 0) {
			mergedAlignments.addAll(existingAlignments);
		}
		mergedAlignments.addAll(newAlignments);
		if (mergedAlignments.size() == 0) {
			builder.rmAttribute(VcfInfoAttributes.BREAKEND_ALIGNMENTS.attribute());
		} else {
			builder.attribute(VcfInfoAttributes.BREAKEND_ALIGNMENTS.attribute(), mergedAlignments);
		}
		nextRecord = builder.make();
	}
	private List<String> writeAlignmentAnnotation(List<SAMRecord> alignments) {
		List<String> aln = new ArrayList<>(alignments.size());
		for (SAMRecord r : alignments) {
			String sb = r.getReferenceName().replace('|', '_') +
					':' +
					r.getAlignmentStart() +
					'|' +
					(r.getReadNegativeStrandFlag() ? '-' : '+') +
					'|' +
					r.getCigarString() +
					'|' +
					r.getMappingQuality();
			aln.add(sb);
			// TODO: better bwa XA tag validation
			if (r.hasAttribute("XA")) {
				// (chr,STRANDpos,CIGAR,NM;)*
				for (String xs : r.getAttribute("XA").toString().split(";")) {
					String[] fields = xs.split(",");
					if(fields.length == 4 && fields[1].length() > 1) {
						String xa = fields[0] +
								':' +
								fields[1].substring(1) +
								'|' +
								fields[1].charAt(0) +
								'|' +
								fields[2] +
								'|';
						aln.add(xa);
					}
				}
			}
		}
		return aln;
	}
	@Override
	public void close() {
		log.debug("Closing UntemplatedSequenceAnnotator");
		// TODO: close feeding thread more cleanly than just shutting down the process
		vcfStream.close();
		//alignerStream.close();
	}
}
