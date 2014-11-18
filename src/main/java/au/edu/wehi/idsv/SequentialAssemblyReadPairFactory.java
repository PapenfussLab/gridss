package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;

import com.google.common.collect.PeekingIterator;

/**
 * Iterates over GRIDSS assemblies
 * @author Daniel Cameron
 *
 */
public class SequentialAssemblyReadPairFactory extends SequentialSAMRecordFactoryBase<SAMRecord> {
	/**
	 * <p>mate iterator <b>must</b> be mate-coordinate sorted. @see SAMRecordMateCoordinateComparator<p>
	 * @param reads
	 * @param mate
	 */
	public SequentialAssemblyReadPairFactory(
			PeekingIterator<SAMRecord> mates) {
		super(mates);
	}
	public SAMRecordAssemblyEvidence createAssembly(SAMRecord record, ProcessingContext processContext, AssemblyEvidenceSource source) {
		if (record.getReadUnmappedFlag()) return null;
		SAMRecord mate = findAssociatedSAMRecord(record);
		SAMRecord assembly;
		SAMRecord realign;
		boolean recordIsAssembly = record.getFirstOfPairFlag(); 
		if (recordIsAssembly) {
			assembly = record;
			realign = mate;
		} else {
			assembly = record;
			realign = mate;
		}
		SAMRecordAssemblyEvidence evidence = AssemblyFactory.incorporateRealignment(processContext, new SAMRecordAssemblyEvidence(source, assembly, realign), realign);
		if (!recordIsAssembly) {
			if (evidence instanceof RealignedSAMRecordAssemblyEvidence) {
				evidence = new RealignedRemoteSAMRecordAssemblyEvidence(processContext, source, mate, record);
			} else {
				// realignment not good enough to qualify as breakend
				evidence = null;
			}
		}
		return evidence;
	}
	@Override
	public SAMRecord findAssociatedSAMRecord(SAMRecord record) {
		return findMatching(record.getReferenceIndex(), record.getAlignmentStart(), record.getReadName() + (record.getFirstOfPairFlag() ? "/1" : "/2"));
	}
	@Override
	protected String getMatchingKey(SAMRecord record) {
		return record.getReadName() + (record.getFirstOfPairFlag() ? "/2" : "/1");
	}
	@Override
	protected int getReferenceIndex(SAMRecord record) {
		return record.getMateReferenceIndex();
	}
	@Override
	protected int getPosition(SAMRecord record) {
		return record.getMateAlignmentStart();
	}
}
