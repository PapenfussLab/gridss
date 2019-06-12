package scambler;

import au.edu.wehi.idsv.LinearGenomicCoordinate;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.visualisation.GexfHelper;
import com.google.common.collect.Maps;
import com.google.common.collect.Ordering;
import com.google.common.primitives.Doubles;
import htsjdk.samtools.SAMRecord;
import it.uniroma1.dis.wsngroup.gexf4j.core.*;
import it.uniroma1.dis.wsngroup.gexf4j.core.data.Attribute;
import it.uniroma1.dis.wsngroup.gexf4j.core.data.AttributeClass;
import it.uniroma1.dis.wsngroup.gexf4j.core.data.AttributeList;
import it.uniroma1.dis.wsngroup.gexf4j.core.data.AttributeType;
import it.uniroma1.dis.wsngroup.gexf4j.core.impl.GexfImpl;
import it.uniroma1.dis.wsngroup.gexf4j.core.impl.data.AttributeListImpl;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

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
	public void exportOverlapGraph(List<SAMRecord> reads, int minOverlap, File file, LinearGenomicCoordinate lgc) {
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
		Attribute attrEliminated = edgeAttrList.createAttribute("eliminated", AttributeType.BOOLEAN, "eliminated");
		// Attribute attrOverlap = edgeAttrList.createAttribute("l", AttributeType.INTEGER, "overlap length");
		// Attribute attrMatch = edgeAttrList.createAttribute("m", AttributeType.INTEGER, "base matches");
		// Attribute attrMismatch = edgeAttrList.createAttribute("mm", AttributeType.INTEGER, "base mismatches");
		graph.getAttributeLists().add(edgeAttrList);
		
		HashMap<Read, Node> startLookup = Maps.newHashMap();
		HashMap<Read, Node> endLookup = Maps.newHashMap();
		OverlapLookup ol = new OverlapLookup(minOverlap);
		for (SAMRecord sam : reads) {
			Read r = Read.create(lgc, sam);
			Node startnode = graph.createNode("start_" + sam.getReadName() + "/" + (SAMRecordUtil.getSegmentIndex(sam) + 1));
			Node endnode = graph.createNode("end_" + sam.getReadName() + "/" + (SAMRecordUtil.getSegmentIndex(sam) + 1));
			startLookup.put(r, startnode);
			endLookup.put(r, endnode);
			ol.add(r);
			r.sanityCheck();
		}
		for (Read r : startLookup.keySet()) {
			Node rstart = startLookup.get(r);
			Node rend = endLookup.get(r);
			for (Overlap o : ol.successors(r)) {
				Node ostart = startLookup.get(o.read2);
				Node oend = endLookup.get(o.read2);
				
				Edge commonEdge = ostart.connectTo(rend).setEdgeType(EdgeType.DIRECTED);
				commonEdge.setWeight(o.overlap);
				commonEdge.getAttributeValues().createValue(attrSeq, new String(o.read2.getSeq().getBytes(0, o.overlap)));
				
				Edge beforeEdge = rstart.connectTo(ostart).setEdgeType(EdgeType.DIRECTED);
				beforeEdge.setWeight(r.getSeq().length() - o.overlap);
				beforeEdge.getAttributeValues().createValue(attrSeq, new String(o.read1.getSeq().getBytes(0, o.read1.getSeq().length() - o.overlap)));
				
				Edge afterEdge = rend.connectTo(oend).setEdgeType(EdgeType.DIRECTED);
				afterEdge.setWeight(o.read2.getSeq().length() - o.overlap);
				afterEdge.getAttributeValues().createValue(attrSeq, new String(o.read2.getSeq().getBytes(o.overlap, o.read2.getSeq().length() - o.overlap)));
			}
		}
		transitiveReduce(graph, attrEliminated);
		GexfHelper.saveTo(gexf, file);
	}
	private static final Ordering<Edge> byEdgeWeight= new Ordering<Edge>() {
		@Override
		public int compare(Edge left, Edge right) {
			return Doubles.compare(left.getWeight(), right.getWeight());
		}
	};
	private void transitiveReduce(Graph graph, Attribute attrEliminated) {
		// TODO: pre-sort the out edges so we don't need to do it every time
		for (Node v : graph.getNodes()) {
			List<Edge> vedges = v.getEdges().stream()
				.filter(e -> e.getSource() == v)
				.sorted(byEdgeWeight)
				.collect(Collectors.toList());
			Set<Node> inplay = vedges.stream()
				.map(e -> e.getTarget())
				.collect(Collectors.toSet());
			Set<Node> eliminated = new HashSet<>();
			for (Edge vw : vedges) {
				Node w = vw.getTarget();
				if (inplay.contains(w)) {
					List<Edge> wedges = w.getEdges().stream()
							.filter(e -> e.getSource() == w)
							.sorted(byEdgeWeight)
							.collect(Collectors.toList());
					for (Edge wx : wedges) {
						Node x = wx.getTarget();
						if (inplay.contains(x)) {
							inplay.remove(x);
							eliminated.add(x);
						}
					}
				}
			}
			// skipping Myers 2005 Fig 4. lines 15-19 for now
			for (Edge vw : vedges) {
				Node w = vw.getTarget();
				if (eliminated.contains(w)) {
					vw.getAttributeValues().createValue(attrEliminated, "true");
				}
			}
		}
	}
	/**
	 * Compresses layout graph
	 */
	public void compress() {
		
	}
}
