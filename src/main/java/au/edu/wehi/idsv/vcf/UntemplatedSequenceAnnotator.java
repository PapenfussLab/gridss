package au.edu.wehi.idsv.vcf;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

import au.edu.wehi.idsv.GenomicProcessingContext;
import au.edu.wehi.idsv.IdsvVariantContext;
import au.edu.wehi.idsv.IdsvVariantContextBuilder;
import au.edu.wehi.idsv.VariantContextDirectedEvidence;
import au.edu.wehi.idsv.alignment.ExternalProcessStreamingAligner;
import au.edu.wehi.idsv.util.AutoClosingIterator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class UntemplatedSequenceAnnotator implements CloseableIterator<IdsvVariantContext> {
	public static final byte DEFAULT_QUAL_SCORE = 20;
	private static final Log log = Log.getInstance(UntemplatedSequenceAnnotator.class);
	private final GenomicProcessingContext context;
	private final File vcf;
	private final List<String> cmd;
	private final int threads;
	private CloseableIterator<IdsvVariantContext> vcfStream;
	private PeekingIterator<SAMRecord> samStream;
	private ExternalProcessStreamingAligner aligner;
	private Thread feedingAligner;
	private IdsvVariantContext nextRecord = null;
	public UntemplatedSequenceAnnotator(GenomicProcessingContext context, File vcf, List<String> aligner_command_line, int threads) {
		this.context = context;
		this.vcf = vcf;
		this.cmd = aligner_command_line;
		this.threads = threads;
		this.vcfStream = getVcf();
		createRecordAlignmentStream();
		this.samStream = Iterators.peekingIterator(this.aligner);
	}
	private CloseableIterator<IdsvVariantContext> getVcf() {
		VCFFileReader vcfReader = new VCFFileReader(vcf, false);
		CloseableIterator<VariantContext> it = vcfReader.iterator();
		Iterator<IdsvVariantContext> idsvIt = Iterators.transform(it, variant -> IdsvVariantContext.create(context, null, variant));
		return new AutoClosingIterator<IdsvVariantContext>(idsvIt, it, vcfReader);
	}
	private ExternalProcessStreamingAligner createRecordAlignmentStream() {
		aligner = new ExternalProcessStreamingAligner(context.getSamReaderFactory(), cmd, context.getReferenceFile(), threads);
		feedingAligner = new Thread(() -> feedExternalAligner());
		feedingAligner.setName("usa-to-aligner");
		feedingAligner.start();
		return aligner;
	}
	private void feedExternalAligner() {
		try (CloseableIterator<IdsvVariantContext> it = getVcf()) {
			while (it.hasNext()) {
				IdsvVariantContext vc = it.next();
				if (vc instanceof VariantContextDirectedEvidence) {
					VariantContextDirectedEvidence e = (VariantContextDirectedEvidence)vc;
					byte[] seq = e.getBreakendSequence();
					byte[] qual = e.getBreakendQuality();
					if (seq != null && seq.length > 0) {
						if (qual == null) {
							qual = new byte[seq.length];
							Arrays.fill(qual, DEFAULT_QUAL_SCORE);
						}
						FastqRecord fq = new FastqRecord(e.getID(), seq, null, qual);
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
	}
	@Override
	public boolean hasNext() {
		ensureNext();
		return nextRecord != null;
	}
	@Override
	public IdsvVariantContext next() {
		ensureNext();
		if (!hasNext()) throw new NoSuchElementException();
		IdsvVariantContext result = nextRecord;
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
			IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(context, nextRecord);
			builder.attribute(VcfInfoAttributes.BREAKEND_ALIGNMENTS.attribute(), writeAlignmentAnnotation(alignments));
			nextRecord = builder.make();
		}
	}
	private List<String> writeAlignmentAnnotation(List<SAMRecord> alignments) {
		List<String> aln = new ArrayList<>(alignments.size());
		for (SAMRecord r : alignments) {
			StringBuilder sb = new StringBuilder();
			sb.append(r.getReferenceName());
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
