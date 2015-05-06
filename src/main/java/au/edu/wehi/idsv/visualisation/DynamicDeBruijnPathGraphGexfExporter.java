package au.edu.wehi.idsv.visualisation;

import it.uniroma1.dis.wsngroup.gexf4j.core.Edge;
import it.uniroma1.dis.wsngroup.gexf4j.core.EdgeType;
import it.uniroma1.dis.wsngroup.gexf4j.core.Gexf;
import it.uniroma1.dis.wsngroup.gexf4j.core.Graph;
import it.uniroma1.dis.wsngroup.gexf4j.core.IDType;
import it.uniroma1.dis.wsngroup.gexf4j.core.IntervalType;
import it.uniroma1.dis.wsngroup.gexf4j.core.Mode;
import it.uniroma1.dis.wsngroup.gexf4j.core.Node;
import it.uniroma1.dis.wsngroup.gexf4j.core.data.Attribute;
import it.uniroma1.dis.wsngroup.gexf4j.core.data.AttributeClass;
import it.uniroma1.dis.wsngroup.gexf4j.core.data.AttributeList;
import it.uniroma1.dis.wsngroup.gexf4j.core.data.AttributeType;
import it.uniroma1.dis.wsngroup.gexf4j.core.dynamic.Spell;
import it.uniroma1.dis.wsngroup.gexf4j.core.dynamic.TimeFormat;
import it.uniroma1.dis.wsngroup.gexf4j.core.impl.GexfImpl;
import it.uniroma1.dis.wsngroup.gexf4j.core.impl.SpellImpl;
import it.uniroma1.dis.wsngroup.gexf4j.core.impl.data.AttributeListImpl;

import java.io.File;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import au.edu.wehi.idsv.debruijn.DeBruijnPathGraph;
import au.edu.wehi.idsv.debruijn.DeBruijnPathNode;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;

public class DynamicDeBruijnPathGraphGexfExporter<T, PN extends DeBruijnPathNode<T>> implements DeBruijnPathGraphExporter<T, PN> {
	//private static final Log log = Log.getInstance(DeBruijnPathGraphGexfExporter.class);
	private final HashMap<PN, Node> lookup = Maps.newHashMap();
	private final Gexf gexf = new GexfImpl();
	private final Graph graph;
	private final Attribute attrTotalWeight;
	private final Attribute attrMaxKmerWeight;
	private final Attribute attrLength;
	private final Attribute attrIsReference;
	private final Attribute attrContig;
	private final Attribute attrContigList;
	private final Attribute attrSubgraphId;
	private final Attribute attrStartNode;
	private int currentTime = 0;
	public DynamicDeBruijnPathGraphGexfExporter() {
		gexf.getMetadata()
			.setLastModified(new Date())
			.setCreator("GRIDSS")
			.setDescription("de Bruijn kmer path graph. Non-branching adjacent kmer graph nodes have been merged.");
		gexf.setVisualization(true);
		graph = gexf.getGraph()
			.setIDType(IDType.STRING)
			.setTimeType(TimeFormat.INTEGER)
			.setDefaultEdgeType(EdgeType.DIRECTED)
			.setMode(Mode.DYNAMIC);
			
		AttributeList nodeAttrList = new AttributeListImpl(AttributeClass.NODE)
			.setMode(Mode.DYNAMIC);
		attrTotalWeight = nodeAttrList.createAttribute("tw", AttributeType.INTEGER, "total weight");
		attrMaxKmerWeight = nodeAttrList.createAttribute("mw", AttributeType.INTEGER, "max kmer weight");
		graph.getAttributeLists().add(nodeAttrList);
		
		AttributeList nodeStaticAttrList = new AttributeListImpl(AttributeClass.NODE)
			.setMode(Mode.STATIC);
		attrLength = nodeStaticAttrList.createAttribute("l", AttributeType.INTEGER, "length of path");
		attrIsReference = nodeStaticAttrList.createAttribute("r", AttributeType.BOOLEAN, "isReference");
		attrSubgraphId = nodeStaticAttrList.createAttribute("sg", AttributeType.INTEGER, "subgraph index");
		attrStartNode = nodeStaticAttrList.createAttribute("sn", AttributeType.BOOLEAN, "starting node");
		attrContigList = nodeStaticAttrList.createAttribute("contigs", AttributeType.LISTSTRING, "Contig memberships");
		attrContig = nodeStaticAttrList.createAttribute("contig", AttributeType.BOOLEAN, "Part of assembly");
		graph.getAttributeLists().add(nodeStaticAttrList);
		
		//AttributeList edgeAttrList = new AttributeListImpl(AttributeClass.EDGE)
		//	.setMode(Mode.STATIC);
		//edgeAttrList.createAttribute("b", AttributeType.STRING, "base");
		//graph.getAttributeLists().add(edgeAttrList);
	}
	@Override
	public DeBruijnPathGraphExporter<T, PN> snapshot(DeBruijnPathGraph<T, PN> pg) {
		currentTime++;
		for (PN pn : pg.getPaths()) {
			if (!lookup.containsKey(pn)) {
				// Add node
				Node node = graph.createNode(pn.toString(pg.getGraph()));
				node.getAttributeValues().createValue(attrLength, ((Integer)pn.getPath().size()).toString());
				node.getAttributeValues().createValue(attrIsReference, ((Boolean)pg.isReference(pn)).toString());
				lookup.put(pn, node);
			}
			trackChanges(pg, pn);
			
		}
		for (PN pn : pg.getPaths()) {
			ensureNextEdges(pg, pn);
		}
		HashSet<PN> allNodes = new HashSet<PN>(pg.allNodes());
		// set node validity timing
		for (PN pn : lookup.keySet()) {
			ensureNodeState(lookup.get(pn), allNodes.contains(pn));
		}
		return this;
	}
	private void ensureNodeState(Node node, boolean active) {
		if (!active) {
			if (node.getSpells().size() > 0) {
				Spell lastSpell = node.getSpells().get(node.getSpells().size() - 1);
				if (!lastSpell.hasEndDate()) {
					lastSpell.setEndValue(currentTime);
					lastSpell.setEndIntervalType(IntervalType.OPEN);
				}
			}
		} else {
			if (node.getSpells().size() > 0 && !node.getSpells().get(node.getSpells().size() - 1).hasEndDate()) {
				// already active - nothing to do
			} else {
				Spell s = new SpellImpl();
				s.setStartValue(currentTime);
				s.setStartIntervalType(IntervalType.CLOSE);
				node.getSpells().add(s);
			}
		}
	}
	private void trackChanges(DeBruijnPathGraph<T, PN> pg, PN pn) {
		Node node = lookup.get(pn);
		GexfHelper.setDynamicAttribute(node, currentTime, attrTotalWeight, pn.weight());
		GexfHelper.setDynamicAttribute(node, currentTime, attrMaxKmerWeight, pn.maxNodeWeight(pg.getGraph()));
	}
	private void ensureNextEdges(DeBruijnPathGraph<T, PN> pg, PN pn) {
		Node node = lookup.get(pn);
		for (PN next : pg.next(pn)) {
			Node nextNode = lookup.get(next);
			ensureEdge(node, nextNode);
		}
	}
	private Edge ensureEdge(Node from, Node to) {
		for (Edge e : from.getEdges()) {
			if (e.getTarget() == to) return e;
		}
		return from.connectTo(to).setEdgeType(EdgeType.DIRECTED);
	}
	/* (non-Javadoc)
	 * @see au.edu.wehi.idsv.visualisation.DeBruijnPathGraphExporter#contig(java.util.List)
	 */
	@Override
	public DeBruijnPathGraphExporter<T, PN> contigs(List<List<PN>> assembledContigs) {
		currentTime++;
		Multimap<PN, Integer> mm = ArrayListMultimap.create();
		int contig = 0;
		for (List<PN> assembledContig : assembledContigs) {
			for (PN pn : assembledContig) {
				mm.put(pn, contig);
			}
		}
		for (PN pn : mm.keySet()) {
			Node node = lookup.get(pn);
			node.getAttributeValues().createValue(attrContig, "true");
			StringBuilder sb = new StringBuilder();
			for (Integer i : mm.get(pn)) {
				sb.append(';');
				sb.append(i);
			}
			sb.replace(0, 1, ""); // strip starting seperator
			node.getAttributeValues().createValue(attrContigList, sb.toString());
		}
		return this;
	}
	/* (non-Javadoc)
	 * @see au.edu.wehi.idsv.visualisation.DeBruijnPathGraphExporter#saveTo(java.io.File)
	 */
	@Override
	public DeBruijnPathGraphExporter<T, PN> saveTo(File file) {
		GexfHelper.saveTo(gexf, file);
		return this;
	}
	@Override
	public void annotateSubgraphs(List<Set<PN>> subgraphs) {
		currentTime++;
		int i = 0;
		for (Set<PN> sg : subgraphs) {
			for (PN pn : sg) {
				lookup.get(pn).getAttributeValues()
					.createValue(attrSubgraphId, ((Integer)i).toString());
			}
			i++;
		}
	}
	@Override
	public void annotateStartingPaths(List<Iterable<PN>> startingPaths) {
		currentTime++;
		for (Iterable<PN> sg : startingPaths) {
			for (PN pn : sg) {
				lookup.get(pn).getAttributeValues()
					.createValue(attrStartNode, "true");
			}
		}
	}
}
