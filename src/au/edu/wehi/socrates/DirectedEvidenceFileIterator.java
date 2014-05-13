package au.edu.wehi.socrates;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;

import htsjdk.samtools.util.Log;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.util.CloseableIterator;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import com.google.common.collect.Iterators;

public class DirectedEvidenceFileIterator implements CloseableIterator<DirectedEvidence> {
	private final SAMFileReader svReader;
	private final SAMFileReader mateReader;
	private final SAMFileReader realignReader;
	private final VCFFileReader vcfReader;
	private final SAMRecordIterator svIt;
	private final SAMRecordIterator mateIt;
	private final SAMRecordIterator realignIt;
	private final CloseableIterator<VariantContext> vcfIt;
	private final DirectedEvidenceIterator it;
	private final Log log = Log.getInstance(DirectedEvidenceFileIterator.class);
	public DirectedEvidenceFileIterator(
			ProcessingContext processContext,
			File sv,
			File mate,
			File realign,
			File vcf) {
		svReader = sv == null ? null : new SAMFileReader(sv);
		mateReader = mate == null ? null : new SAMFileReader(mate);
		realignReader = realign == null ? null : new SAMFileReader(sv);
		vcfReader = vcf == null ? null : new VCFFileReader(vcf);;
		svIt = svReader == null ? null : svReader.iterator();
		mateIt = mateReader == null ? null : mateReader.iterator();
		realignIt = realignReader == null ? null : realignReader.iterator();
		vcfIt = vcfReader == null ? null : vcfReader.iterator();
		it = new DirectedEvidenceIterator(
				processContext,
				Iterators.peekingIterator(svIt),
				Iterators.peekingIterator(mateIt),
				Iterators.peekingIterator(realignIt),
				Iterators.peekingIterator(vcfIt));
	}
	@Override
	public void close() {
		close(svIt, mateIt, realignIt, vcfIt);
		close(svReader, mateReader, realignReader, vcfReader);
	}
	private void close(Closeable... toClose) {
		for (Closeable c : toClose) {
			if (c != null) {
				try {
					c.close();
				} catch (IOException e) {
					// log and swallow
					log.warn(e);
				}
			}
		}
	}
	@Override
	public boolean hasNext() {
		return it.hasNext();
	}
	@Override
	public DirectedEvidence next() {
		return it.next();
	}
	@Deprecated
	@Override
	public void remove() {
		it.remove();
	}
}
