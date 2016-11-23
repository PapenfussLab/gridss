package au.edu.wehi.idsv.visualisation;

import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.FileFilter;
import java.io.IOException;

import org.apache.commons.io.filefilter.WildcardFileFilter;
import org.junit.Ignore;
import org.junit.Test;
import org.junit.experimental.categories.Category;

import com.google.common.collect.ImmutableList;

import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.Hg19Tests;
import au.edu.wehi.idsv.IntermediateFilesTest;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.SAMEvidenceSource;
import au.edu.wehi.idsv.sam.SAMFileUtil;
import au.edu.wehi.idsv.util.FixMissingMateInformation;
import gridss.ComputeSamTags;
import htsjdk.samtools.SAMFileHeader.SortOrder;


public class StaticDeBruijnPathGraphGexfExporterTest extends IntermediateFilesTest {
	@Test
	@Category(Hg19Tests.class)
	public void positional_should_export_dot() throws IOException {
		File output = new File(super.testFolder.getRoot(), "chr12-244000.vcf");
		setReference(Hg19Tests.findHg19Reference("chr12.fa"));
		createInput(new File("src/test/resources/chr12-244000.tagged.bam"));
		ProcessingContext pc = getCommandlineContext();
		SAMEvidenceSource ses = new SAMEvidenceSource(pc, input, 0);
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, ImmutableList.of(ses), output);
		aes.assembleBreakends(null);
		File dir = new File(super.testFolder.getRoot(), "visualisation");
		File[] export = dir.listFiles((FileFilter)new WildcardFileFilter("*.dot"));
		// TODO: change graph write timing
		assertTrue(export != null && export.length > 0);
	}
	@Override
	public ProcessingContext getCommandlineContext() {
		ProcessingContext pc = super.getCommandlineContext();
		pc.getAssemblyParameters().errorCorrection.maxBaseMismatchForCollapse = 1;
		pc.getAssemblyParameters().errorCorrection.collapseBubblesOnly = true;
		pc.getConfig().getVisualisation().assemblyGraph = true;
		pc.getConfig().getVisualisation().directory.mkdirs();
		return pc;
	}
	@Test
	//@Ignore("Once-off data import")
	public void add_tags() throws IOException {
		File in = new File("src/test/resources/chr12-244000.bam");
		File insq = new File("src/test/resources/chr12-244000.sq.bam");
		File outsq = new File("src/test/resources/chr12-244000.sq.tagged.bam");
		File outfixmate = new File("src/test/resources/chr12-244000.sq.tagged.fixmate.bam");
		File out = new File("src/test/resources/chr12-244000.tagged.bam");
		SAMFileUtil.sort(getFSContext(), in, insq, SortOrder.queryname);
		File ref = Hg19Tests.findHg19Reference("chr12.fa");
		ComputeSamTags cmd = new ComputeSamTags();
		cmd.instanceMain(new String[] {
				"INPUT=" + insq.getAbsolutePath(),
				"OUTPUT=" + outsq.getAbsolutePath(),
				"REFERENCE_SEQUENCE=" + ref.getAbsolutePath()
		});
		FixMissingMateInformation cmd2 = new FixMissingMateInformation();
		cmd2.instanceMain(new String[] {
				"INPUT=" + outsq.getAbsolutePath(),
				"OUTPUT=" + outfixmate.getAbsolutePath()
		});
		SAMFileUtil.sort(getFSContext(), outfixmate, out, SortOrder.coordinate);
	}
}
