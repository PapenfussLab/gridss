package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import au.edu.wehi.idsv.sam.SAMRecordUtil;

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
			assembly = mate;
			realign = record;
		}
		SAMRecordAssemblyEvidence evidence = AssemblyFactory.incorporateRealignment(processContext, new SAMRecordAssemblyEvidence(source, assembly, realign), realign);
		if (!recordIsAssembly) {
			if (evidence instanceof RealignedSAMRecordAssemblyEvidence) {
				evidence = ((RealignedSAMRecordAssemblyEvidence)evidence).asRemote();
			} else {
				// realignment not good enough to qualify as breakend
				evidence = null;
			}
		}
		return evidence;
	}
	@Override
	public SAMRecord findAssociatedSAMRecord(SAMRecord record) {
		if (record == null) return null;
		if (record.getReadUnmappedFlag()) return null;
		return findMatching(record.getReferenceIndex(), record.getAlignmentStart(), record.getReadName() + (record.getFirstOfPairFlag() ? SAMRecordUtil.FIRST_OF_PAIR_NAME_SUFFIX : SAMRecordUtil.SECOND_OF_PAIR_NAME_SUFFIX));
	}
	@Override
	protected String getMatchingKey(SAMRecord record) {
		return record.getReadName() + (record.getFirstOfPairFlag() ? SAMRecordUtil.SECOND_OF_PAIR_NAME_SUFFIX : SAMRecordUtil.FIRST_OF_PAIR_NAME_SUFFIX);
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
