package au.edu.wehi.idsv.vcf;

import au.edu.wehi.idsv.*;
import au.edu.wehi.idsv.alignment.ExternalProcessStreamingAligner;
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
	private final File referenceFile;
	private final File vcf;
	private final boolean overwrite;
	private final List<String> cmd;
	private final int threads;
	private final SAMSequenceDictionary dict;
	private CloseableIterator<VariantContext> vcfStream;
	private PeekingIterator<SAMRecord> samStream;
	private ExternalProcessStreamingAligner aligner;
	private Thread feedingAligner;
	private VariantContext nextRecord = null;
	public UntemplatedSequenceAnnotator(File referenceFile, File vcf, boolean overwrite, List<String> aligner_command_line, int threads, SAMSequenceDictionary dict) {
		this.referenceFile = referenceFile;
		this.vcf = vcf;
		this.overwrite = overwrite;
		this.cmd = aligner_command_line;
		this.threads = threads;
		this.vcfStream = getVcf();
		createRecordAlignmentStream();
		this.samStream = Iterators.peekingIterator(this.aligner);
		this.dict = dict;
	}
	private CloseableIterator<VariantContext> getVcf() {
		VCFFileReader vcfReader = new VCFFileReader(vcf, false);
		CloseableIterator<VariantContext> it = vcfReader.iterator();
		return new AutoClosingIterator<VariantContext>(it, vcfReader);
	}
	private ExternalProcessStreamingAligner createRecordAlignmentStream() {
		log.debug("Started async external alignment feeder thread.");
		aligner = new ExternalProcessStreamingAligner(SamReaderFactory.make(), cmd, referenceFile, threads, dict);
		feedingAligner = new Thread(() -> feedExternalAligner());
		feedingAligner.setName("async-feedExternalAligner");
		feedingAligner.start();
		return aligner;
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

	private void feedExternalAligner() {
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
			if (aligner.hasNext()) {
				log.debug("Traversing aligner stream to enable graceful termination.");
				while(aligner.hasNext()) {
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
		while (samStream.hasNext() && samStream.peek().getReadName().equals(readName)) {
			SAMRecord r = samStream.next();
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
		try {
			aligner.close();
		} catch (IOException e) {
			log.warn(e);
		}
	}
}
