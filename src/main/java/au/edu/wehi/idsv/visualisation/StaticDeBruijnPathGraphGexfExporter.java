package au.edu.wehi.idsv.visualisation;

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

import java.io.File;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Set;

import au.edu.wehi.idsv.debruijn.DeBruijnNodeBase;
import au.edu.wehi.idsv.debruijn.DeBruijnPathGraph;
import au.edu.wehi.idsv.debruijn.PathNode;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;

public class StaticDeBruijnPathGraphGexfExporter<T extends DeBruijnNodeBase, PN extends PathNode<T>>
	implements DeBruijnPathGraphExporter<T, PN> {
	//private static final Log log = Log.getInstance(DeBruijnPathGraphGexfExporter.class);
	private final int k;
	private final HashMap<PN, Node> lookup = Maps.newHashMap();
	private final Gexf gexf = new GexfImpl();
	private final Graph graph;
	private final Attribute attrTotalWeight;
	private final Attribute attrMaxKmerWeight;
	private final Attribute attrLength;
	private final Attribute attrContig;
	private final Attribute attrContigList;
	private final Attribute attrSubgraphId;
	private final Attribute attrStartNode;
	private final AttributeList nodeStaticAttrList;
	public StaticDeBruijnPathGraphGexfExporter(int k) {
		this.k = k;
		gexf.getMetadata()
			.setLastModified(new Date())
			.setCreator("GRIDSS")
			.setDescription("de Bruijn kmer path graph");
		gexf.setVisualization(true);
		graph = gexf.getGraph()
			.setIDType(IDType.STRING)
			.setDefaultEdgeType(EdgeType.DIRECTED)
			.setMode(Mode.STATIC);
		nodeStaticAttrList = new AttributeListImpl(AttributeClass.NODE)
			.setMode(Mode.STATIC);
		attrLength = nodeStaticAttrList.createAttribute("l", AttributeType.INTEGER, "length of path");
		attrSubgraphId = nodeStaticAttrList.createAttribute("sg", AttributeType.INTEGER, "subgraph index");
		attrStartNode = nodeStaticAttrList.createAttribute("sn", AttributeType.BOOLEAN, "starting node");
		attrContigList = nodeStaticAttrList.createAttribute("contigs", AttributeType.LISTSTRING, "Contig memberships");
		attrContig = nodeStaticAttrList.createAttribute("contig", AttributeType.BOOLEAN, "Part of assembly");
		attrTotalWeight = nodeStaticAttrList.createAttribute("tw", AttributeType.INTEGER, "total weight");
		attrMaxKmerWeight = nodeStaticAttrList.createAttribute("mw", AttributeType.INTEGER, "max kmer weight");
		graph.getAttributeLists().add(nodeStaticAttrList);
	}
	protected AttributeList getNodeStaticAttributeList() { return nodeStaticAttrList; }
	@Override
	public DeBruijnPathGraphExporter<T, PN> snapshot(DeBruijnPathGraph<T, PN> pg) {
		if (lookup.size() != 0) throw new IllegalStateException("Cannot add more than one snapshot");
		for (PN pn : pg.getPaths()) {
			if (!lookup.containsKey(pn)) {
				// Add node
				Node node = graph.createNode(pn.pathKmerString(k));
				setStaticAttributes(node, pn);
				lookup.put(pn, node);
			}
		}
		for (PN pn : pg.getPaths()) {
			ensureNextEdges(pg, pn);
		}
		return this;
	}
	protected void setStaticAttributes(Node node, PN pn) {
		node.getAttributeValues().createValue(attrLength, ((Integer)pn.getPath().size()).toString());
		node.getAttributeValues().createValue(attrTotalWeight, ((Integer)pn.getWeight()).toString());
		node.getAttributeValues().createValue(attrMaxKmerWeight, ((Integer)pn.getMaxKmerWeight()).toString());
	}
	private void ensureNextEdges(DeBruijnPathGraph<T, PN> pg, PN pn) {
		Node node = lookup.get(pn);
		for (PN next : pg.nextPath(pn)) {
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
	public void annotateStartingPaths(List<Set<PN>> startingPaths) {
		for (Set<PN> sg : startingPaths) {
			for (PN pn : sg) {
				lookup.get(pn).getAttributeValues()
					.createValue(attrStartNode, "true");
			}
		}
	}
}
