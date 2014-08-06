package au.edu.wehi.idsv;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;

import java.io.File;
import java.util.Iterator;
import java.util.List;

import au.edu.wehi.idsv.sam.SAMFileUtil;

import com.google.common.collect.Lists;

public class SortRealignedSoftClips {
	private static final Log log = Log.getInstance(SortRealignedSoftClips.class);
	private final ProcessingContext processContext;
	private final SAMEvidenceSource source;
	private final List<SAMFileWriter> scwriters = Lists.newArrayList();
	private final List<SAMFileWriter> realignmentWriters = Lists.newArrayList();
	public SortRealignedSoftClips(ProcessingContext processContext, SAMEvidenceSource source) {
		this.processContext = processContext;
		this.source = source;
	}
	public boolean isComplete() {
		FileSystemContext fsc = processContext.getFileSystemContext();
		boolean done = true;
		if (processContext.shouldProcessPerChromosome()) {
			for (SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
				done &= fsc.getSoftClipRemoteBamForChr(source.getSourceFile(), seq.getSequenceName()).exists();
				done &= fsc.getRealignmentRemoteBamForChr(source.getSourceFile(), seq.getSequenceName()).exists();
			}
		} else {
			done &= fsc.getSoftClipRemoteBam(source.getSourceFile()).exists();
			done &= fsc.getRealignmentRemoteBam(source.getSourceFile()).exists();
		}
		return done;
	}
	public void process(boolean force) {
		if (force) {
			deleteOutput();
		}
		if (isComplete()) {
			log.debug("SortRealignedSoftClips: no work to do for ", source.getSourceFile());
			return;
		}
		if (!source.isRealignmentComplete()) {
			String msg = String.format("Soft clip realignment for %s not completed. Unable to process", source.getSourceFile());
			log.error(msg);
			throw new IllegalStateException(msg);
		}
		try {
			log.info("START: sorting mapped soft clips for ", source.getSourceFile());
			createUnsortedOutputWriters();
			writeUnsortedOutput();
			closeUnsortedWriters();
			sort();
			deleteUnsorted();
			log.info("SUCCESS: sorting mapped soft clips for ", source.getSourceFile());
		} catch (Exception e) {
			String msg = String.format("Unable to sort mapped soft clips for %s", source.getSourceFile());
			log.error(e, msg);
			closeUnsortedWriters();
			deleteOutput();
			throw new RuntimeException(msg, e);
		}
	}
	private void deleteOutput() {
		deleteUnsorted();
		FileSystemContext fsc = processContext.getFileSystemContext();
		for (SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
			tryDelete(fsc.getSoftClipRemoteBamForChr(source.getSourceFile(), seq.getSequenceName()));
			tryDelete(fsc.getRealignmentRemoteBamForChr(source.getSourceFile(), seq.getSequenceName()));
		}
		tryDelete(fsc.getSoftClipRemoteUnsortedBam(source.getSourceFile()));
		tryDelete(fsc.getRealignmentRemoteUnsortedBam(source.getSourceFile()));
	}
	private void deleteUnsorted() {
		FileSystemContext fsc = processContext.getFileSystemContext();
		for (SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
			tryDelete(fsc.getSoftClipRemoteUnsortedBamForChr(source.getSourceFile(), seq.getSequenceName()));
			tryDelete(fsc.getRealignmentRemoteUnsortedBamForChr(source.getSourceFile(), seq.getSequenceName()));
		}
		tryDelete(fsc.getSoftClipRemoteUnsortedBam(source.getSourceFile()));
		tryDelete(fsc.getRealignmentRemoteUnsortedBam(source.getSourceFile()));
	}
	private void tryDelete(File f) {
		try {
			if (f.exists()) {
				if (!f.delete()) {
					log.error("Unable to delete intermediate file ", f,  " during rollback.");
				}
			}
		} catch (Exception e) {
			log.error(e, "Unable to delete intermediate file ", f,  " during rollback.");
		}
	}
	private void sort() {
		FileSystemContext fsc = processContext.getFileSystemContext();
		if (processContext.shouldProcessPerChromosome()) {
			for (SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
				SAMFileUtil.sort(processContext,
						fsc.getSoftClipRemoteUnsortedBamForChr(source.getSourceFile(), seq.getSequenceName()),
						fsc.getSoftClipRemoteBamForChr(source.getSourceFile(), seq.getSequenceName()),
						new RealignedSoftClipEvidence.RealignmentCoordinateComparator());
				SAMFileUtil.sort(processContext,
						fsc.getRealignmentRemoteUnsortedBamForChr(source.getSourceFile(), seq.getSequenceName()),
						fsc.getRealignmentRemoteBamForChr(source.getSourceFile(), seq.getSequenceName()),
						SortOrder.coordinate);
			}
		} else {
			SAMFileUtil.sort(processContext,
					fsc.getSoftClipRemoteUnsortedBam(source.getSourceFile()),
					fsc.getSoftClipRemoteBam(source.getSourceFile()),
					new RealignedSoftClipEvidence.RealignmentCoordinateComparator());
			SAMFileUtil.sort(processContext,
					fsc.getRealignmentRemoteUnsortedBam(source.getSourceFile()),
					fsc.getRealignmentRemoteBam(source.getSourceFile()),
					SortOrder.coordinate);
		}
	}
	private void closeUnsortedWriters() {
		for (SAMFileWriter w : scwriters) {
			CloserUtil.close(w);
		}
		scwriters.clear();
		for (SAMFileWriter w : realignmentWriters) {
			CloserUtil.close(w);
		}
		realignmentWriters.clear();
	}
	private void writeUnsortedOutput() {
		Iterator<DirectedEvidence> it = source.iterator();
		while (it.hasNext()) {
			DirectedEvidence de = it.next();
			if (de instanceof RealignedSoftClipEvidence) {
				RealignedSoftClipEvidence evidence = (RealignedSoftClipEvidence)de;
				scwriters.get(evidence.getBreakendSummary().referenceIndex2 % scwriters.size()).addAlignment(evidence.getSAMRecord());
				realignmentWriters.get(evidence.getBreakendSummary().referenceIndex2 % realignmentWriters.size()).addAlignment(evidence.getRealignedSAMRecord());
			}
		}
	}
	private void createUnsortedOutputWriters() {
		FileSystemContext fsc = processContext.getFileSystemContext();
		SAMFileHeader header = source.getHeader().clone();
		header.setSortOrder(SortOrder.unsorted);
		if (processContext.shouldProcessPerChromosome()) {
			for (SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
				scwriters.add(processContext.getSamReaderWriterFactory().makeSAMOrBAMWriter(header, true, fsc.getSoftClipRemoteUnsortedBamForChr(source.getSourceFile(), seq.getSequenceName())));
				realignmentWriters.add(processContext.getSamReaderWriterFactory().makeSAMOrBAMWriter(header, true, fsc.getRealignmentRemoteUnsortedBamForChr(source.getSourceFile(), seq.getSequenceName())));
			}
		} else {
			scwriters.add(processContext.getSamReaderWriterFactory().makeSAMOrBAMWriter(header, true, fsc.getSoftClipRemoteUnsortedBam(source.getSourceFile())));
			realignmentWriters.add(processContext.getSamReaderWriterFactory().makeSAMOrBAMWriter(header, true, fsc.getRealignmentRemoteUnsortedBam(source.getSourceFile())));
		}
	}
}
