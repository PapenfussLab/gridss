package au.edu.wehi.idsv;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.PeekingIterator;

import au.edu.wehi.idsv.sam.SAMRecordUtil;
import htsjdk.samtools.SAMRecord;

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
		SAMRecord mate = findFirstAssociatedSAMRecord(record);
		SAMRecord assembly;
		SAMRecord realign;
		boolean isLocalAssembly = record.getFirstOfPairFlag(); 
		if (isLocalAssembly) {
			assembly = record;
			realign = mate;
		} else {
			assembly = mate;
			realign = record;
		}
		SAMRecordAssemblyEvidence evidence = AssemblyFactory.hydrate(source, assembly);
		if (realign != null) {
			evidence = AssemblyFactory.incorporateRealignment(processContext, evidence, ImmutableList.of(realign));
		}
		if (!isLocalAssembly && evidence instanceof RealignedSAMRecordAssemblyEvidence) {
			evidence = ((RealignedSAMRecordAssemblyEvidence)evidence).asRemote();
		}
		return evidence;
	}
	public SAMRecord findFirstAssociatedSAMRecord(SAMRecord record) {
		if (record == null) return null;
		if (record.getReadUnmappedFlag()) return null;
		return findFirstMatching(record.getReferenceIndex(), record.getAlignmentStart(), record.getReadName() + (record.getFirstOfPairFlag() ? SAMRecordUtil.FIRST_OF_PAIR_NAME_SUFFIX : SAMRecordUtil.SECOND_OF_PAIR_NAME_SUFFIX));
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
	@Override
	protected String trackedName() {
		return "matepairing.arp";
	}
}
