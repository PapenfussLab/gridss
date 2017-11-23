package scambler;

import java.io.File;
import java.util.Date;
import java.util.HashMap;
import java.util.List;

import com.google.common.collect.Maps;

import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.visualisation.GexfHelper;
import htsjdk.samtools.SAMRecord;
import it.uniroma1.dis.wsngroup.gexf4j.core.Edge;
import it.uniroma1.dis.wsngroup.gexf4j.core.EdgeType;
import it.uniroma1.dis.wsngroup.gexf4j.core.Gexf;
import it.uniroma1.dis.wsngroup.gexf4j.core.Graph;
import it.uniroma1.dis.wsngroup.gexf4j.core.IDType;
import it.uniroma1.dis.wsngroup.gexf4j.core.Mode;
import it.uniroma1.dis.wsngroup.gexf4j.core.Node;
import it.uniroma1.dis.wsngroup.gexf4j.core.data.Attribute;
import it.uniroma1.dis.wsngroup.gexf4j.core.data.AttributeClass;
import it.uniroma1.dis.wsngroup.gexf4j.core.data.AttributeList;
import it.uniroma1.dis.wsngroup.gexf4j.core.data.AttributeType;
import it.uniroma1.dis.wsngroup.gexf4j.core.impl.GexfImpl;
import it.uniroma1.dis.wsngroup.gexf4j.core.impl.data.AttributeListImpl;

/**
 * Proof of concept static graph implementation
 * 
 * @author Daniel Cameron
 *
 */
public class Poc {
	// OverlapGraph { read, read, offset }
	// ReadPairGraph { read1, read2 }
	// FragmentSizeDistribution
	
	// CompressedOverlapNode { _seq, (read, contig_offset) }
	
	// Algorithm details to resolve:
	// - how to allocate read overlap
	// 		- just allocate reads at contig resolution
	// - How to generate haplotypes?
	// - efficient overlap calculation
	// - how to deal with errors
	// 		- sequencing errors
	//		- chimeric fragments
	// Am I just making an OLC haplotype caller? (would that be problematic anyway?)
	
	/**
	 * Load overlap string graph
	 */
	public void exportOverlapGraph(List<SAMRecord> reads, int minOverlap, File file) {
		Gexf gexf = new GexfImpl();
			gexf.getMetadata()
			.setLastModified(new Date())
			.setCreator("GRIDSS")
			.setDescription("Raw overlap string graph");
		gexf.setVisualization(true);
		Graph graph = gexf.getGraph()
			.setIDType(IDType.STRING)
			.setDefaultEdgeType(EdgeType.DIRECTED)
			.setMode(Mode.STATIC);
		AttributeList nodeAttrList = new AttributeListImpl(AttributeClass.NODE)
			.setMode(Mode.STATIC);
		graph.getAttributeLists().add(nodeAttrList);
		AttributeList edgeAttrList = new AttributeListImpl(AttributeClass.EDGE)
			.setMode(Mode.STATIC);
		Attribute attrSeq = edgeAttrList.createAttribute("seq", AttributeType.STRING, "sequence");
		// Attribute attrOverlap = edgeAttrList.createAttribute("l", AttributeType.INTEGER, "overlap length");
		// Attribute attrMatch = edgeAttrList.createAttribute("m", AttributeType.INTEGER, "base matches");
		// Attribute attrMismatch = edgeAttrList.createAttribute("mm", AttributeType.INTEGER, "base mismatches");
		graph.getAttributeLists().add(edgeAttrList);
		
		HashMap<Read, Node> startLookup = Maps.newHashMap();
		HashMap<Read, Node> endLookup = Maps.newHashMap();
		OverlapLookup ol = new OverlapLookup(minOverlap);
		for (SAMRecord sam : reads) {
			Read r = Read.create(sam);
			Node startnode = graph.createNode("start_" + sam.getReadName() + "/" + (SAMRecordUtil.getSegmentIndex(sam) + 1));
			Node endnode = graph.createNode("end_" + sam.getReadName() + "/" + (SAMRecordUtil.getSegmentIndex(sam) + 1));
			startLookup.put(r, startnode);
			endLookup.put(r, endnode);
			ol.add(r);
			r.sanityCheck();
		}
		for (Read r : startLookup.keySet()) {
			Node rend = endLookup.get(r);
			for (Overlap o : ol.successors(r)) {
				Node ostart = startLookup.get(o.read2);
				Edge edge = ostart.connectTo(rend).setEdgeType(EdgeType.DIRECTED);
				edge.setWeight(o.overlap);
				edge.getAttributeValues().createValue(attrSeq, new String(o.read2.getSeq().getBytes(0, o.overlap)));
			}
		}
		GexfHelper.saveTo(gexf, file);
	}
	/**
	 * Compresses layout graph
	 */
	public void compress() {
		
	}
}
