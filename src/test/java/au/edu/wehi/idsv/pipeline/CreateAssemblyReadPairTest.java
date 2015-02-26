package au.edu.wehi.idsv.pipeline;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.SequenceUtil;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;

import org.junit.Before;
import org.junit.Test;

import au.edu.wehi.idsv.AssemblyEvidence;
import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.AssemblyFactory;
import au.edu.wehi.idsv.BreakendSummary;
import au.edu.wehi.idsv.BreakpointFastqEncoding;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.IntermediateFilesTest;
import au.edu.wehi.idsv.MockReadEvidenceAssembler;
import au.edu.wehi.idsv.ProcessStep;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.ReadEvidenceAssembler;
import au.edu.wehi.idsv.RealignedSAMRecordAssemblyEvidence;
import au.edu.wehi.idsv.RemoteEvidence;
import au.edu.wehi.idsv.SAMEvidenceSource;
import au.edu.wehi.idsv.SAMRecordAssemblyEvidence;
import au.edu.wehi.idsv.SmallIndelSAMRecordAssemblyEvidence;
import au.edu.wehi.idsv.SmallIndelSAMRecordAssemblyEvidenceTest;
import au.edu.wehi.idsv.util.AutoClosingIterator;

import com.google.common.base.Predicate;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;


public class CreateAssemblyReadPairTest extends IntermediateFilesTest {
	@Before
	public void setup() throws IOException {
		super.setup();
		evidence = Lists.newArrayList();
		realign = Lists.newArrayList();
	}
	List<SAMRecordAssemblyEvidence> evidence;
	List<SAMRecord> realign;
	int realignCount;
	private void orderedAddNoRealign(SAMRecordAssemblyEvidence e) {
		SAMRecord r = new SAMRecord(getContext().getBasicSamHeader());
		r.setReadUnmappedFlag(true);
		orderedAdd(e, r);
	}
	private void orderedAdd(SAMRecordAssemblyEvidence e, int referenceIndex, int position, boolean negativeStrand) {
		SAMRecord r = new SAMRecord(getContext().getBasicSamHeader());
		r.setReferenceIndex(referenceIndex);
		r.setAlignmentStart(position);
		r.setReadUnmappedFlag(false);
		r.setCigar(new Cigar(ImmutableList.of(new CigarElement(e.getBreakendSequence().length, CigarOperator.MATCH_OR_MISMATCH))));
		r.setMappingQuality(30);
		orderedAdd(e, r);
	}
	private void orderedAdd(SAMRecordAssemblyEvidence e, SAMRecord r) {
		evidence.add(e);
		if (r != null) {
			r.setReadName(BreakpointFastqEncoding.getRealignmentFastq(e).getReadHeader());
			byte[] seq = e.getBreakendSequence();
			if (r.getReadNegativeStrandFlag()) {
				SequenceUtil.reverseComplement(seq);
			}
			r.setReadBases(e.getBreakendSequence());
		}
		realign.add(r);
	}
	public static class CompletedSAMEvidenceSource extends MockSAMEvidenceSource {
		public CompletedSAMEvidenceSource() {
			super(IntermediateFilesTest.getContext(), 0, 300);
		}
		@Override
		public boolean isComplete(ProcessStep step) {
			return true;
		}
		@Override
		protected CloseableIterator<DirectedEvidence> iterator(
				final boolean includeReadPair,
				final boolean includeSoftClip,
				final boolean includeSoftClipRemote,
				final File readPair,
				final File pairMate,
				final File softClip,
				final File realigned,
				final File remoteSoftClip,
				final File remoteRealigned,
				final String chr) {
			return new AutoClosingIterator<DirectedEvidence>(Collections.<DirectedEvidence>emptyIterator());
		}
	}
	/**
	 * Dirty hack instead of proper dependency injection.
	 * @author Daniel Cameron
	 *
	 */
	public static class FixedAssemblyEvidenceSource extends AssemblyEvidenceSource {
		public FixedAssemblyEvidenceSource(ProcessingContext processContext, List<SAMRecordAssemblyEvidence> assemblyResult, File file) {
			super(processContext, ImmutableList.<SAMEvidenceSource>of(new CompletedSAMEvidenceSource()), file);
			this.assemblies = Lists.newArrayList();
			if (processContext.shouldProcessPerChromosome()) {
				// will fail when multithreaded
				for (final SAMSequenceRecord chr : processContext.getDictionary().getSequences()) {
					assemblies.add(Lists.newArrayList(Iterables.filter(assemblyResult, new Predicate<SAMRecordAssemblyEvidence>() {
						@Override
						public boolean apply(SAMRecordAssemblyEvidence input) {
							return input.getBreakendSummary().referenceIndex == chr.getSequenceIndex();
						}
					})));
				}
			} else {
				assemblies.add(assemblyResult);
			}
		}
		private List<List<SAMRecordAssemblyEvidence>> assemblies;
		int calls = 0;
		@SuppressWarnings("unchecked")
		@Override
		protected ReadEvidenceAssembler getAssembler() {
			// works for one call if not per chr, otherwise returns results for chr in order
			return new MockReadEvidenceAssembler((List<AssemblyEvidence>)(Object)assemblies.get(calls++));
		}
		
	}
	private void writeRealign(ProcessingContext pc, File base) {
		realignCount = 0;
		if (pc.shouldProcessPerChromosome()) {
			for (SAMSequenceRecord chr : SMALL_FA.getSequenceDictionary().getSequences()) {
				SAMFileWriter writer = pc.getSamFileWriterFactory(false).makeSAMOrBAMWriter(pc.getBasicSamHeader(), true,
						pc.getFileSystemContext().getRealignmentBamForChr(base, chr.getSequenceName()));
				for (int i = 0; i < evidence.size(); i++) {
					if (realign.get(i) != null && evidence.get(i).getBreakendSummary().referenceIndex == chr.getSequenceIndex()) {
						writer.addAlignment(realign.get(i));
						if (!realign.get(i).getReadUnmappedFlag()) realignCount++;
					}
				}
				writer.close();
			}
		} else {
			SAMFileWriter writer = pc.getSamFileWriterFactory(false).makeSAMOrBAMWriter(pc.getBasicSamHeader(), true,
					pc.getFileSystemContext().getRealignmentBam(base));
			for (int i = 0; i < evidence.size(); i++) {
				if (realign.get(i) != null) {
					writer.addAlignment(realign.get(i));
					if (!realign.get(i).getReadUnmappedFlag()) realignCount++;
				}
			}
			writer.close();
		}
		for (SAMRecordAssemblyEvidence e : evidence) {
			if (e instanceof SmallIndelSAMRecordAssemblyEvidence) {
				// considered realigned already
				realignCount++;
			}
		}
	}
	private void go() {
		go(getCommandlineContext(false));
		go(getCommandlineContext(true));
	}
	private void go(ProcessingContext pc) { go(pc, true); }
	private void go(ProcessingContext pc, boolean writefiltered) {
		pc.getAssemblyParameters().writeFilteredAssemblies = writefiltered;
		pc.getAssemblyParameters().minReads = 0;
		pc.getRealignmentParameters().minLength = 0;
		pc.getRealignmentParameters().mapqUniqueThreshold = 0;
		pc.getRealignmentParameters().minAverageQual = 0;
		File faes = new File(output + "-realign-" + pc.shouldProcessPerChromosome());
		File frp = new File(output + "-readpair-" + pc.shouldProcessPerChromosome());
		writeAssemblies(pc, faes);
		writeAssemblies(pc, frp);
		writeRealign(pc, faes);
		writeRealign(pc, frp);
		AssemblyEvidenceSource ar = new AssemblyEvidenceSource(pc, ImmutableList.<SAMEvidenceSource>of(SES()), faes);
		AssemblyEvidenceSource rp = new AssemblyEvidenceSource(pc, ImmutableList.<SAMEvidenceSource>of(SES()), frp);
		rp.ensureAssembled();
		assertMatch(rp, ar);
	}
	private void writeAssemblies(ProcessingContext pc, File file) {
		new FixedAssemblyEvidenceSource(pc, evidence, file).ensureAssembled();
	}
	private void assertMatch(AssemblyEvidenceSource rp, AssemblyEvidenceSource aes) {
		assertEquals(evidence.size(), Iterators.size(aes.iterator(false, true)));
		assertEquals(evidence.size(), Iterators.size(rp.iterator(false, true)));
		assertEquals(evidence.size() + realignCount, Iterators.size(rp.iterator(true, true)));
		assertContains(Lists.newArrayList(rp.iterator(false, true)), Lists.newArrayList(aes.iterator(false, true)));
		assertSorted(Lists.newArrayList(rp.iterator(false, true)));
		assertSorted(Lists.newArrayList(aes.iterator(false, true)));
		for (SAMSequenceRecord chr : SMALL_FA.getSequenceDictionary().getSequences()) {
			assertContains(Lists.newArrayList(rp.iterator(false, true, chr.getSequenceName())), Lists.newArrayList(aes.iterator(false, true, chr.getSequenceName())));
			assertCorrectChr(Lists.newArrayList(rp.iterator(false, true, chr.getSequenceName())), chr.getSequenceIndex());
			assertCorrectChr(Lists.newArrayList(rp.iterator(false, true, chr.getSequenceName())), chr.getSequenceIndex());
			assertSorted(Lists.newArrayList(rp.iterator(false, true, chr.getSequenceName())));
			assertSorted(Lists.newArrayList(aes.iterator(false, true, chr.getSequenceName())));
		}
		// Remote realignments exist
		for (RealignedSAMRecordAssemblyEvidence e : Lists.newArrayList(Iterators.filter(rp.iterator(true, true), RealignedSAMRecordAssemblyEvidence.class))) {
			BreakendSummary remote = e.getBreakendSummary().remoteBreakpoint();
			assertContains(Lists.newArrayList(rp.iterator(true, true)), remote);
			assertContains(Lists.newArrayList(rp.iterator(true, true, SMALL_FA.getSequenceDictionary().getSequences().get(remote.referenceIndex).getSequenceName())), remote);
			if (e instanceof RemoteEvidence) {
				assertContains(Lists.newArrayList(aes.iterator(false, true)), e.getBreakendSummary().remoteBreakpoint());
			}
		}
	}
	private void assertCorrectChr(Iterable<SAMRecordAssemblyEvidence> list, int expected) {
		for (SAMRecordAssemblyEvidence s : list) {
			assertEquals(expected, s.getBreakendSummary().referenceIndex);
		}
	}
	private void assertSorted(Iterable<SAMRecordAssemblyEvidence> list) {
		BreakendSummary last = new BreakendSummary(0,  FWD, 1,  1);
		for (SAMRecordAssemblyEvidence s : list) {
			assertTrue(BreakendSummary.ByStartEnd.compare(last, s.getBreakendSummary()) <= 0);
			last = s.getBreakendSummary();
		}
	}
	private void assertContains(Iterable<SAMRecordAssemblyEvidence> superset, Iterable<SAMRecordAssemblyEvidence> subset) {
		for (SAMRecordAssemblyEvidence s : subset) {
			assertContains(superset, s);
		}
	}
	private void assertContains(Iterable<SAMRecordAssemblyEvidence> list, SAMRecordAssemblyEvidence r) {
		for (SAMRecordAssemblyEvidence s : list) {
			if (s.getEvidenceID().equals(r.getEvidenceID()) && s.getBreakendSummary().equals(r.getBreakendSummary())) return;
		}
		fail(String.format("Could not find %s", r.getEvidenceID()));
	}
	private void assertContains(Iterable<SAMRecordAssemblyEvidence> list, BreakendSummary bs) {
		for (SAMRecordAssemblyEvidence s : list) {
			if (bs.equals(s.getBreakendSummary())) return;
		}
		fail(String.format("Could not find %s", bs));
	}
	@Test
	public void read_pair_should_match_assembly_realign_iterator() {
		ProcessingContext pc = getContext();
		AssemblyEvidenceSource aes = AES(pc);
		orderedAddNoRealign(AssemblyFactory.createAnchored(pc, aes, FWD, Sets.<DirectedEvidence>newHashSet(), 0, 1, 1, B("AA"), B("AA"), 0, 0));
		orderedAddNoRealign(AssemblyFactory.createAnchored(pc, aes, BWD, Sets.<DirectedEvidence>newHashSet(), 0, 1, 1, B("AA"), B("AA"), 0, 0));
		orderedAddNoRealign(AssemblyFactory.createAnchored(pc, aes, FWD, Sets.<DirectedEvidence>newHashSet(), 0, 2, 1, B("AA"), B("AA"), 0, 0));
		orderedAddNoRealign(AssemblyFactory.createAnchored(pc, aes, BWD, Sets.<DirectedEvidence>newHashSet(), 0, 2, 1, B("AA"), B("AA"), 0, 0));
		orderedAdd(AssemblyFactory.createAnchored(pc, aes, BWD, Sets.<DirectedEvidence>newHashSet(), 0, 3, 1, B("AA"), B("AA"), 0, 0), 1, 2, true);
		orderedAdd(AssemblyFactory.createAnchored(pc, aes, BWD, Sets.<DirectedEvidence>newHashSet(), 0, 4, 1, B("AA"), B("AA"), 0, 0), 1, 2, false);
		orderedAdd(AssemblyFactory.createAnchored(pc, aes, FWD, Sets.<DirectedEvidence>newHashSet(), 0, 5, 1, B("AA"), B("AA"), 0, 0), 1, 2, true);
		orderedAdd(AssemblyFactory.createAnchored(pc, aes, FWD, Sets.<DirectedEvidence>newHashSet(), 0, 6, 1, B("AA"), B("AA"), 0, 0), 1, 2, false);
		orderedAddNoRealign(AssemblyFactory.createAnchored(pc, aes, FWD, Sets.<DirectedEvidence>newHashSet(), 1, 1, 1, B("AA"), B("AA"), 0, 0));
		orderedAddNoRealign(AssemblyFactory.createAnchored(pc, aes, BWD, Sets.<DirectedEvidence>newHashSet(), 1, 1, 1, B("AA"), B("AA"), 0, 0));
		orderedAddNoRealign(AssemblyFactory.createAnchored(pc, aes, FWD, Sets.<DirectedEvidence>newHashSet(), 1, 2, 1, B("AA"), B("AA"), 0, 0));
		orderedAddNoRealign(AssemblyFactory.createAnchored(pc, aes, BWD, Sets.<DirectedEvidence>newHashSet(), 1, 2, 1, B("AA"), B("AA"), 0, 0));
		orderedAdd(AssemblyFactory.createAnchored(pc, aes, BWD, Sets.<DirectedEvidence>newHashSet(), 1, 3, 1, B("AA"), B("AA"), 0, 0), 1, 2, true);
		orderedAdd(AssemblyFactory.createAnchored(pc, aes, BWD, Sets.<DirectedEvidence>newHashSet(), 1, 4, 1, B("AA"), B("AA"), 0, 0), 1, 2, false);
		orderedAdd(AssemblyFactory.createAnchored(pc, aes, FWD, Sets.<DirectedEvidence>newHashSet(), 1, 5, 1, B("AA"), B("AA"), 0, 0), 1, 2, true);
		orderedAdd(AssemblyFactory.createAnchored(pc, aes, FWD, Sets.<DirectedEvidence>newHashSet(), 1, 6, 1, B("AA"), B("AA"), 0, 0), 1, 2, false);
		go();
	}
	@Test
	public void should_have_realign() {
		orderedAdd(AssemblyFactory.createAnchored(getContext(), AES(), FWD, Sets.<DirectedEvidence>newHashSet(),
				0, 6, 1, B("TT"), B("TT"), 0, 0), 1, 2, false);
		go();
	}
	@Test
	public void read_pair_should_match_assembly_realign_iterator_simple() {
		orderedAddNoRealign(AssemblyFactory.createAnchored(getContext(), AES(), FWD, Sets.<DirectedEvidence>newHashSet(), 0, 1, 1, B("AA"), B("AA"), 0, 0));
		go();
	}
	@Test
	public void assembly_evidence_source_should_resort_to_evidence_order() {
		orderedAddNoRealign(AssemblyFactory.createAnchored(getContext(), AES(), FWD, Sets.<DirectedEvidence>newHashSet(), 0, 10, 5, B("AAAAAAAAAA"), B("AAAAAAAAAA"), 0, 0)); // alignment starts at 5
		orderedAddNoRealign(AssemblyFactory.createAnchored(getContext(), AES(), BWD, Sets.<DirectedEvidence>newHashSet(), 0, 6, 5, B("AAAAAAAAAA"), B("AAAAAAAAAA"), 0, 0)); // alignment starts at 6
		go();
	}
	@Test
	public void perChr_should_write_single_file_for_debugging_purposes() {
		go();
		assertTrue(getCommandlineContext(true).getFileSystemContext().getAssembly(new File(output + "-readpair-true")).exists());
	}
	@Test
	public void should_filter_breakpoints() {
		orderedAdd(AssemblyFactory.createAnchored(getContext(), AES(), FWD, Sets.<DirectedEvidence>newHashSet(),
				0, 5, 1, B("TT"), B("TT"), 0, 0),
				0, 6, false);
		
		ProcessingContext pc = getCommandlineContext(false);
		pc.getAssemblyParameters().writeFilteredAssemblies = false;
		pc.getAssemblyParameters().minReads = 0;
		pc.getRealignmentParameters().minLength = 0;
		pc.getRealignmentParameters().mapqUniqueThreshold = 0;
		pc.getRealignmentParameters().minAverageQual = 0;
		File faes = new File(output + "-realign-" + pc.shouldProcessPerChromosome());
		File frp = new File(output + "-readpair-" + pc.shouldProcessPerChromosome());
		writeAssemblies(pc, faes);
		writeAssemblies(pc, frp);
		writeRealign(pc, faes);
		writeRealign(pc, frp);
		AssemblyEvidenceSource ar = new AssemblyEvidenceSource(pc, ImmutableList.<SAMEvidenceSource>of(SES()), faes);
		AssemblyEvidenceSource rp = new AssemblyEvidenceSource(pc, ImmutableList.<SAMEvidenceSource>of(SES()), frp);
		rp.ensureAssembled();
		assertEquals("precondition: breakend and realign should have been written as breakend passes filters", 1, Iterators.size(ar.iterator(false,  true)));
		List<SAMRecordAssemblyEvidence> result = Lists.newArrayList(rp.iterator(false,  false));
		assertEquals(0, result.size());
	}
	@Test
	public void should_include_small_indel_remote_breakends() {
		orderedAddNoRealign(SmallIndelSAMRecordAssemblyEvidenceTest.create(1, "10M10D10M", "NNNNNNNNNNATATATATAT", FWD));
		go();
	}
}
