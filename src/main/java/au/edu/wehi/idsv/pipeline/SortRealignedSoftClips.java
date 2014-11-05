package au.edu.wehi.idsv.pipeline;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Log;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.Iterator;
import java.util.List;

import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.FileSystemContext;
import au.edu.wehi.idsv.ProcessStep;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.RealignedRemoteSoftClipEvidence;
import au.edu.wehi.idsv.RealignedSoftClipEvidence;
import au.edu.wehi.idsv.SAMEvidenceSource;
import au.edu.wehi.idsv.sam.SAMFileUtil;
import au.edu.wehi.idsv.sam.SAMFileUtil.SortCallable;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;

public class SortRealignedSoftClips extends DataTransformStep {
	private static final Log log = Log.getInstance(SortRealignedSoftClips.class);
	private final SAMEvidenceSource source;
	private final List<SAMFileWriter> scwriters = new ArrayList<>();
	private final List<SAMFileWriter> realignmentWriters = new ArrayList<>();
	public SortRealignedSoftClips(ProcessingContext processContext, SAMEvidenceSource source) {
		super(processContext);
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
	@Override
	public void process(EnumSet<ProcessStep> steps) {
		if (isComplete()) {
			log.debug("SortRealignedSoftClips: no work to do for ", source.getSourceFile());
		}
		if (!canProcess()) {
			String msg = String.format("Soft clip realignment for %s not completed. Unable to process", source.getSourceFile());
			log.error(msg);
			throw new IllegalStateException(msg);
			// return EnumSet.of(ProcessStep.SORT_REALIGNED_SOFT_CLIPS);
		}
		try {
			log.info("START: sorting mapped soft clips for ", source.getSourceFile());
			createUnsortedOutputWriters();
			writeUnsortedOutput();
			close();
			sort();
			deleteTemp();
			log.info("SUCCESS: sorting mapped soft clips for ", source.getSourceFile());
		} catch (Exception e) {
			String msg = String.format("Unable to sort mapped soft clips for %s", source.getSourceFile());
			log.error(e, msg);
			close();
			deleteTemp();
			deleteOutput();
			throw new RuntimeException(msg, e);
			// return EnumSet.of(ProcessStep.SORT_REALIGNED_SOFT_CLIPS);
		}
	}
	private void sort() throws IOException {
		FileSystemContext fsc = processContext.getFileSystemContext();
		List<SortCallable> actions = new ArrayList<>();
		if (processContext.shouldProcessPerChromosome()) {
			for (SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
				actions.add(new SAMFileUtil.SortCallable(processContext,
						fsc.getSoftClipRemoteUnsortedBamForChr(source.getSourceFile(), seq.getSequenceName()),
						fsc.getSoftClipRemoteBamForChr(source.getSourceFile(), seq.getSequenceName()),
						new RealignedSoftClipEvidence.RealignmentCoordinateComparator()));
				actions.add(new SAMFileUtil.SortCallable(processContext,
						fsc.getRealignmentRemoteUnsortedBamForChr(source.getSourceFile(), seq.getSequenceName()),
						fsc.getRealignmentRemoteBamForChr(source.getSourceFile(), seq.getSequenceName()),
						SortOrder.coordinate));
			}
		} else {
			actions.add(new SAMFileUtil.SortCallable(processContext,
					fsc.getSoftClipRemoteUnsortedBam(source.getSourceFile()),
					fsc.getSoftClipRemoteBam(source.getSourceFile()),
					new RealignedSoftClipEvidence.RealignmentCoordinateComparator()));
			actions.add(new SAMFileUtil.SortCallable(processContext,
					fsc.getRealignmentRemoteUnsortedBam(source.getSourceFile()),
					fsc.getRealignmentRemoteBam(source.getSourceFile()),
					SortOrder.coordinate));
		}
		// PARALLEL opportunity - not great candidate due memory usage of sorting
		for (SortCallable c : actions) {
			c.call();
		}
	}
	@Override
	public void close() {
		super.close();
		scwriters.clear();
		realignmentWriters.clear();
	}
	private void writeUnsortedOutput() {
		Iterator<RealignedSoftClipEvidence> it = Iterators.filter(source.iterator(), RealignedSoftClipEvidence.class);
		while (it.hasNext()) {
			DirectedEvidence de = it.next();
			if (de instanceof RealignedRemoteSoftClipEvidence) continue;
			RealignedSoftClipEvidence evidence = (RealignedSoftClipEvidence)de;
			scwriters.get(evidence.getBreakendSummary().referenceIndex2 % scwriters.size()).addAlignment(evidence.getSAMRecord());
			realignmentWriters.get(evidence.getBreakendSummary().referenceIndex2 % realignmentWriters.size()).addAlignment(evidence.getRealignedSAMRecord());
		}
	}
	private void createUnsortedOutputWriters() {
		FileSystemContext fsc = processContext.getFileSystemContext();
		SAMFileHeader header = source.getHeader().clone();
		header.setSortOrder(SortOrder.unsorted);
		if (processContext.shouldProcessPerChromosome()) {
			for (SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
				scwriters.add(processContext.getSamReaderWriterFactory().makeSAMOrBAMWriter(header, true, fsc.getSoftClipRemoteUnsortedBamForChr(source.getSourceFile(), seq.getSequenceName())));
				toClose.add(scwriters.get(scwriters.size() - 1));
				realignmentWriters.add(processContext.getSamReaderWriterFactory().makeSAMOrBAMWriter(header, true, fsc.getRealignmentRemoteUnsortedBamForChr(source.getSourceFile(), seq.getSequenceName())));
				toClose.add(realignmentWriters.get(realignmentWriters.size() - 1));
			}
		} else {
			scwriters.add(processContext.getSamReaderWriterFactory().makeSAMOrBAMWriter(header, true, fsc.getSoftClipRemoteUnsortedBam(source.getSourceFile())));
			toClose.add(scwriters.get(scwriters.size() - 1));
			realignmentWriters.add(processContext.getSamReaderWriterFactory().makeSAMOrBAMWriter(header, true, fsc.getRealignmentRemoteUnsortedBam(source.getSourceFile())));
			toClose.add(realignmentWriters.get(realignmentWriters.size() - 1));
		}
	}
	@Override
	protected Log getLog() {
		return log;
	}
	@Override
	public List<File> getInputs() {
		return ImmutableList.of();
	}
	@Override
	public boolean canProcess() {
		return source.isComplete(ProcessStep.REALIGN_SOFT_CLIPS);
	}
	@Override
	public List<File> getOutput() {
		List<File> outputs = new ArrayList<>();
		FileSystemContext fsc = processContext.getFileSystemContext();
		if (processContext.shouldProcessPerChromosome()) {
			for (SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
				outputs.add(fsc.getSoftClipRemoteBamForChr(source.getSourceFile(), seq.getSequenceName()));
				outputs.add(fsc.getRealignmentRemoteBamForChr(source.getSourceFile(), seq.getSequenceName()));
			}
		} else {
			outputs.add(fsc.getSoftClipRemoteBam(source.getSourceFile()));
			outputs.add(fsc.getRealignmentRemoteBam(source.getSourceFile()));
		}
		return outputs;
	}
	@Override
	public List<File> getTemporary() {
		List<File> files = new ArrayList<>();
		FileSystemContext fsc = processContext.getFileSystemContext();
		if (processContext.shouldProcessPerChromosome()) {
			for (SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
				files.add(fsc.getSoftClipRemoteUnsortedBamForChr(source.getSourceFile(), seq.getSequenceName()));
				files.add(fsc.getRealignmentRemoteUnsortedBamForChr(source.getSourceFile(), seq.getSequenceName()));
			}
		} else {
			files.add(fsc.getSoftClipRemoteUnsortedBam(source.getSourceFile()));
			files.add(fsc.getRealignmentRemoteUnsortedBam(source.getSourceFile()));
		}
		return files;
	}
}
