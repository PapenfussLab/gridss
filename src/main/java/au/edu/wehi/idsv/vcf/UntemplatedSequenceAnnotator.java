package au.edu.wehi.idsv.vcf;

import au.edu.wehi.idsv.*;
import au.edu.wehi.idsv.alignment.BwaStreamingAligner;
import au.edu.wehi.idsv.alignment.ExternalProcessStreamingAligner;
import au.edu.wehi.idsv.alignment.StreamingAligner;
import au.edu.wehi.idsv.alignment.StreamingAlignerIterator;
import au.edu.wehi.idsv.util.AutoClosingIterator;
import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import joptsimple.internal.Strings;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class UntemplatedSequenceAnnotator implements CloseableIterator<VariantContext> {
	public static final byte DEFAULT_QUAL_SCORE = 20;
	private static final Log log = Log.getInstance(UntemplatedSequenceAnnotator.class);
	private static final Pattern breakendRegex = Pattern.compile("^(.(?<leftins>.*))?[\\[\\]].*[\\[\\]]((?<rightins>.*).)?$");
	private final File vcf;
	private final boolean overwrite;
	private CloseableIterator<VariantContext> vcfStream;
	private PeekingIterator<SAMRecord> alignerStream;
	private Thread feedingAligner;
	private VariantContext nextRecord = null;
	public UntemplatedSequenceAnnotator(File vcf, StreamingAligner aligner, boolean overwrite) {
		this.vcf = vcf;
		this.overwrite = overwrite;
		this.vcfStream = getVcf();
		this.alignerStream = Iterators.peekingIterator(new StreamingAlignerIterator(aligner));
		this.feedingAligner = new Thread(() -> feedStreamingAligner(aligner));
		this.feedingAligner.setName("feedAligner");
		this.feedingAligner.start();
	}

	private CloseableIterator<VariantContext> getVcf() {
		VCFFileReader vcfReader = new VCFFileReader(vcf, false);
		CloseableIterator<VariantContext> it = vcfReader.iterator();
		return new AutoClosingIterator<VariantContext>(it, vcfReader);
	}
	private boolean shouldAttemptAlignment(VariantContext v) {
		return overwrite || !v.hasAttribute(VcfInfoAttributes.BREAKEND_ALIGNMENTS.attribute());
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
	private void feedStreamingAligner(StreamingAligner aligner) {
		try (CloseableIterator<VariantContext> it = getVcf()) {
			while (it.hasNext()) {
				VariantContext vc = it.next();
				if (shouldAttemptAlignment(vc)) {
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
		nextRecord = vcfStream.next();
		annotateNextRecord();
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
		if (alignments.size() > 0) {
			VariantContextBuilder builder = new VariantContextBuilder(nextRecord);
			builder.attribute(VcfInfoAttributes.BREAKEND_ALIGNMENTS.attribute(), writeAlignmentAnnotation(alignments));
			nextRecord = builder.make();
		}
	}
	private List<String> writeAlignmentAnnotation(List<SAMRecord> alignments) {
		List<String> aln = new ArrayList<>(alignments.size());
		for (SAMRecord r : alignments) {
			StringBuilder sb = new StringBuilder();
			sb.append(r.getReferenceName().replace(':', '_').replace('|', '_'));
			sb.append(':');
			sb.append(r.getAlignmentStart());
			sb.append('|');
			sb.append(r.getReadNegativeStrandFlag() ? '-' : '+');
			sb.append('|');
			sb.append(r.getCigarString());
			sb.append('|');
			sb.append(r.getMappingQuality());
			aln.add(sb.toString());
			// TODO: better bwa XA tag validation
			if (r.hasAttribute("XA")) {
				// (chr,STRANDpos,CIGAR,NM;)*
				for (String xs : r.getAttribute("XA").toString().split(";")) {
					String[] fields = xs.split(",");
					if(fields.length == 4 && fields[1].length() > 1) {
						StringBuilder xasb = new StringBuilder();
						xasb.append(fields[0]);
						xasb.append(':');
						xasb.append(fields[1].substring(1));
						xasb.append('|');
						xasb.append(fields[1].charAt(0));
						xasb.append('|');
						xasb.append(fields[2]);
						xasb.append('|');
						aln.add(xasb.toString());
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
