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
import it.uniroma1.dis.wsngroup.gexf4j.core.impl.GexfImpl;
import it.uniroma1.dis.wsngroup.gexf4j.core.impl.SpellImpl;
import it.uniroma1.dis.wsngroup.gexf4j.core.impl.data.AttributeListImpl;

import java.io.File;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Set;

import au.edu.wehi.idsv.debruijn.DeBruijnGraphBase;
import au.edu.wehi.idsv.debruijn.DeBruijnPathGraph;
import au.edu.wehi.idsv.debruijn.subgraph.DeBruijnSubgraphNode;
import au.edu.wehi.idsv.debruijn.subgraph.SubgraphPathNode;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;

public class StaticDeBruijnPathGraphGexfExporter implements DeBruijnPathGraphExporter<DeBruijnSubgraphNode, SubgraphPathNode> {
	//private static final Log log = Log.getInstance(DeBruijnPathGraphGexfExporter.class);
	private final int k;
	private final HashMap<SubgraphPathNode, Node> lookup = Maps.newHashMap();
	private final Gexf gexf = new GexfImpl();
	private final Graph graph;
	private final Attribute attrTotalWeight;
	private final Attribute attrMaxKmerWeight;
	private final Attribute attrLength;
	private final Attribute attrIsReference;
	private final Attribute attrIsNonReference;
	private final Attribute attrContig;
	private final Attribute attrContigList;
	private final Attribute attrSubgraphId;
	private final Attribute attrStartNode;
	private int currentTime = 0;
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
			
		AttributeList nodeStaticAttrList = new AttributeListImpl(AttributeClass.NODE)
			.setMode(Mode.STATIC);
		attrLength = nodeStaticAttrList.createAttribute("l", AttributeType.INTEGER, "length of path");
		attrIsReference = nodeStaticAttrList.createAttribute("r", AttributeType.BOOLEAN, "containsReferenceKmer");
		attrIsNonReference = nodeStaticAttrList.createAttribute("nr", AttributeType.BOOLEAN, "containsNonReferenceKmer");
		attrSubgraphId = nodeStaticAttrList.createAttribute("sg", AttributeType.INTEGER, "subgraph index");
		attrStartNode = nodeStaticAttrList.createAttribute("sn", AttributeType.BOOLEAN, "starting node");
		attrContigList = nodeStaticAttrList.createAttribute("contigs", AttributeType.LISTSTRING, "Contig memberships");
		attrContig = nodeStaticAttrList.createAttribute("contig", AttributeType.BOOLEAN, "Part of assembly");
		attrTotalWeight = nodeStaticAttrList.createAttribute("tw", AttributeType.INTEGER, "total weight");
		attrMaxKmerWeight = nodeStaticAttrList.createAttribute("mw", AttributeType.INTEGER, "max kmer weight");
		graph.getAttributeLists().add(nodeStaticAttrList);
	}
	@Override
	public DeBruijnPathGraphExporter<DeBruijnSubgraphNode, SubgraphPathNode> snapshot(DeBruijnPathGraph<DeBruijnSubgraphNode, SubgraphPathNode> pg) {
		if (lookup.size() != 0) throw new IllegalStateException("Cannot add more than one snapshot");
		currentTime++;
		for (SubgraphPathNode pn : pg.getPaths()) {
			if (!lookup.containsKey(pn)) {
				// Add node
				Node node = graph.createNode(new String(DeBruijnGraphBase.getBaseCalls(pn.getPath(), k)));
				node.getAttributeValues().createValue(attrLength, ((Integer)pn.getPath().size()).toString());
				node.getAttributeValues().createValue(attrIsReference, ((Boolean)pn.containsReferenceKmer()).toString());
				node.getAttributeValues().createValue(attrIsNonReference, ((Boolean)pn.containsNonReferenceKmer()).toString());
				node.getAttributeValues().createValue(attrTotalWeight, ((Integer)pn.getWeight()).toString());
				node.getAttributeValues().createValue(attrMaxKmerWeight, ((Integer)pn.getMaxKmerWeight()).toString());
				lookup.put(pn, node);
			}
		}
		for (SubgraphPathNode pn : pg.getPaths()) {
			ensureNextEdges(pg, pn);
		}
		// set node validity timing
		for (SubgraphPathNode pn : lookup.keySet()) {
			ensureNodeState(lookup.get(pn), pg.getPaths().contains(pn));
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
	private void ensureNextEdges(DeBruijnPathGraph<DeBruijnSubgraphNode, SubgraphPathNode> pg, SubgraphPathNode pn) {
		Node node = lookup.get(pn);
		for (SubgraphPathNode next : pg.nextPath(pn)) {
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
	public DeBruijnPathGraphExporter<DeBruijnSubgraphNode, SubgraphPathNode> contigs(List<List<SubgraphPathNode>> assembledContigs) {
		currentTime++;
		Multimap<SubgraphPathNode, Integer> mm = ArrayListMultimap.create();
		int contig = 0;
		for (List<SubgraphPathNode> assembledContig : assembledContigs) {
			for (SubgraphPathNode pn : assembledContig) {
				mm.put(pn, contig);
			}
		}
		for (SubgraphPathNode pn : mm.keySet()) {
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
	public DeBruijnPathGraphExporter<DeBruijnSubgraphNode, SubgraphPathNode> saveTo(File file) {
		GexfHelper.saveTo(gexf, file);
		return this;
	}
	@Override
	public void annotateSubgraphs(List<Set<SubgraphPathNode>> subgraphs) {
		currentTime++;
		int i = 0;
		for (Set<SubgraphPathNode> sg : subgraphs) {
			for (SubgraphPathNode pn : sg) {
				lookup.get(pn).getAttributeValues()
					.createValue(attrSubgraphId, ((Integer)i).toString());
			}
			i++;
		}
	}
	@Override
	public void annotateStartingPaths(List<Set<SubgraphPathNode>> startingPaths) {
		currentTime++;
		for (Set<SubgraphPathNode> sg : startingPaths) {
			for (SubgraphPathNode pn : sg) {
				lookup.get(pn).getAttributeValues()
					.createValue(attrStartNode, "true");
			}
		}
	}
}
