package assfolder;

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
import scambler.OverlapLookup;

import java.io.File;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

import static java.util.stream.Collectors.groupingBy;

public class OlcExporter {
	public static void exportOverlapGraph(List<SAMRecord> samrecords, int seedLength, int minOverlap, int maxMismatches, File file) {
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
		Attribute attrReadSeq = nodeAttrList.createAttribute("seq", AttributeType.STRING, "sequence");
		graph.getAttributeLists().add(nodeAttrList);
		AttributeList edgeAttrList = new AttributeListImpl(AttributeClass.EDGE)
			.setMode(Mode.STATIC);
		Attribute attrSeq = edgeAttrList.createAttribute("seq", AttributeType.LISTSTRING, "sequence");
		Attribute attrOverlap = edgeAttrList.createAttribute("l", AttributeType.INTEGER, "overlap length");
		Attribute attrOverlapList = edgeAttrList.createAttribute("ll", AttributeType.LISTSTRING, "overlap length");
		Attribute attrMismatch = edgeAttrList.createAttribute("mm", AttributeType.INTEGER, "base mismatches");
		graph.getAttributeLists().add(edgeAttrList);

		OverlapGraph og = new OverlapGraph( seedLength, minOverlap, maxMismatches);
		List<Read> reads = samrecords.stream()
				.map(sr -> new Read(sr))
				.collect(Collectors.toList());
		reads = reads.stream()
				.map(r -> {
					Read dup = og.add(r);
					return r != dup ? null : r;
				})
				.filter(r -> r != null)
				.collect(Collectors.toList());
		reads.stream().forEach(r -> og.add(r));
		Map<String, Node> overlapNodes = reads.stream()
			.map(r -> {
				Node n = graph.createNode(r.uid());
				n.getAttributeValues().createValue(attrReadSeq, new String(r.getBytes(0, r.length())));
				return n;
			})
			.collect(Collectors.toMap(Node::getId, Function.identity()));
		for (Read r : reads) {
			Node n = graph.getNode(r.uid());
			for (List<ReadOverlapSuccessor> list : r.overlapSuccessors.stream().collect(groupingBy(x -> x.read)).values()) {
				list.sort(Comparator.comparing(x -> x.offset));
				ReadOverlapSuccessor ros = list.get(0);
				Edge edge = n.connectTo(graph.getNode(ros.read.uid())).setEdgeType(EdgeType.DIRECTED);
				edge.getAttributeValues().createValue(attrSeq, ros.getLeadingBases() + "/" + ros.getTrailingBases());
				edge.getAttributeValues().createValue(attrOverlap, String.valueOf(ros.overlapLength));
				edge.getAttributeValues().createValue(attrMismatch, String.valueOf(ros.mismatches));
				edge.getAttributeValues().createValue(attrOverlapList, list.stream().map(x -> String.valueOf(x.overlapLength)).collect(Collectors.joining(",")));
			}
		}
		GexfHelper.saveTo(gexf, file);
	}
}
