package au.edu.wehi.idsv.visualisation;

import htsjdk.samtools.util.Log;
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
import it.uniroma1.dis.wsngroup.gexf4j.core.dynamic.TimeFormat;
import it.uniroma1.dis.wsngroup.gexf4j.core.impl.GexfImpl;
import it.uniroma1.dis.wsngroup.gexf4j.core.impl.data.AttributeListImpl;

import java.io.File;
import java.util.Date;
import java.util.HashMap;
import java.util.List;

import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.debruijn.subgraph.DeBruijnSubgraphNode;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class DeBruijnSubgraphGexfExporter implements DeBruijnGraphExporter<DeBruijnSubgraphNode> {
	/**
	 * Maximum size of exported graph.
	 * As soon as the export graph reaches this size, all subsequent graph changes are ignored.
	 */
	private static final int MAXIMUM_NODES = 1000000;
	private static final Log log = Log.getInstance(DeBruijnSubgraphGexfExporter.class);
	private final int k;
	private final HashMap<Long, Node> gexfNodes = Maps.newHashMap();
	private final HashMap<Long, List<Long>> gexfEdges = Maps.newHashMap();
	private final Gexf gexf = new GexfImpl();
	private final Graph graph;
	private final Attribute attrWeight;
	private final Attribute attrIsReference;
	private final Attribute attrReferencePosition;
	private int currentTime;
	public DeBruijnSubgraphGexfExporter(int k) {
		this.k = k;
		gexf.getMetadata()
			.setLastModified(new Date())
			.setCreator("GRIDSS")
			.setDescription("de Bruijn kmer graph");
		gexf.setVisualization(true);
		graph = gexf.getGraph()
			.setIDType(IDType.STRING)
			.setTimeType(TimeFormat.INTEGER)
			.setDefaultEdgeType(EdgeType.DIRECTED)
			.setMode(Mode.DYNAMIC);
			
		AttributeList nodeAttrList = new AttributeListImpl(AttributeClass.NODE)
			.setMode(Mode.DYNAMIC);
		attrWeight = nodeAttrList.createAttribute("w", AttributeType.INTEGER, "weight");
		attrIsReference = nodeAttrList.createAttribute("r", AttributeType.INTEGER, "isReferenceKmer");
		attrReferencePosition = nodeAttrList.createAttribute("rp", AttributeType.INTEGER, "referencePosition");
		graph.getAttributeLists().add(nodeAttrList);
		//AttributeList edgeAttrList = new AttributeListImpl(AttributeClass.EDGE)
		//	.setMode(Mode.STATIC);
		//edgeAttrList.createAttribute("b", AttributeType.STRING, "base");
		//graph.getAttributeLists().add(edgeAttrList);
	}
	@Override
	public DeBruijnSubgraphGexfExporter setTime(int currentTime) {
		this.currentTime = currentTime;
		return this;
	}
	/* (non-Javadoc)
	 * @see au.edu.wehi.idsv.visualisation.DeBruijnGraphExporter#trackChanges(long, au.edu.wehi.idsv.debruijn.DeBruijnNodeBase)
	 */
	@Override
	public DeBruijnSubgraphGexfExporter trackChanges(long kmer, DeBruijnSubgraphNode kmerNode) {
		if (gexfNodes.size() >= MAXIMUM_NODES) return this;
		if (kmerNode == null || kmerNode.getWeight() <= 0) {
			// close out final interval
			Node node = gexfNodes.get(kmer);
			if (node != null) {
				node.getSpells().get(node.getSpells().size() - 1).setEndValue(currentTime);
			}
		} else {
			if (!gexfNodes.containsKey(kmer)) {
				// new node!
				gexfNodes.put(kmer, graph.createNode(KmerEncodingHelper.toString(k, kmer)));
				if (gexfNodes.size() == MAXIMUM_NODES) {
					log.info(String.format("No longer exporting kmer graph updates for this chromosome - max export kmers reached at position %d.", currentTime));
				}
			}
			ensureEdges(kmer);
			Node node = gexfNodes.get(kmer);
			updateDynamicAttributes(node, currentTime, kmerNode);
		}
		return this;
	}
	private void updateDynamicAttributes(Node node, int timeStamp, DeBruijnSubgraphNode kmerNode) {
		GexfHelper.setDynamicAttribute(node, timeStamp, attrWeight, kmerNode.getWeight());
		GexfHelper.setDynamicAttribute(node, timeStamp, attrIsReference, kmerNode.isReference() ? 1 : 0);
		GexfHelper.setDynamicAttribute(node, timeStamp, attrReferencePosition, kmerNode.getExpectedPosition());
		// try getting a sane graph by placing nodes at their anchor position
		//node.setColor(new ColorImpl(kmerNode.isReference() ? 255 : 0, kmerNode.isMateAnchored() ? 255 : 0, 255));
		//node.setPosition(new PositionImpl(
		//		kmerNode.getBestReferencePosition() != null ? kmerNode.getBestReferencePosition() : 
		//			(kmerNode.getMinMatePosition() != null ? kmerNode.getMinMatePosition() : timeStamp), 
		//		kmerNode.isReference() ? 0 : (kmerNode.isMateAnchored() ? 2 : 1),
		//		0));
		//setDynamicAttribute(node, timeStamp, attrMaxReferencePosition, kmerNode.getMaxReferencePosition());
		//setDynamicAttribute(node, timeStamp, attrMinReferencePosition, kmerNode.getMinReferencePosition());
	}
	/**
	 * Ensures that graph edges exist between all adjacent kmers
	 * @param state
	 */
	private void ensureEdges(long state) {
		for (Long next : KmerEncodingHelper.nextStates(k, state)) {
			ensureEdge(state, next);			
		}
		for (Long prev : KmerEncodingHelper.prevStates(k, state)) {
			ensureEdge(prev, state);
		}
	}
	private void ensureEdge(
			long sourceKmer,
			long destKmer) {
		Node sourceNode = gexfNodes.get(sourceKmer);
		if (sourceNode == null) return;
		Node destNode = gexfNodes.get(destKmer);
		if (destNode == null) return;
		List<Long> sourceEdgesAlreadyAdded = gexfEdges.get(sourceKmer);
		if (sourceEdgesAlreadyAdded == null) {
			gexfEdges.put(sourceKmer, Lists.<Long>newArrayList());
			sourceEdgesAlreadyAdded = gexfEdges.get(sourceKmer);
		}
		if (sourceEdgesAlreadyAdded.contains(destKmer)) return;
		Edge edge = sourceNode.connectTo(destNode);
		String kmerBases = KmerEncodingHelper.toString(k, destKmer);
		edge.setLabel(kmerBases.substring(kmerBases.length() - 1));
		edge.setEdgeType(EdgeType.DIRECTED);
		sourceEdgesAlreadyAdded.add(destKmer);
	}
	/* (non-Javadoc)
	 * @see au.edu.wehi.idsv.visualisation.DeBruijnGraphExporter#saveTo(java.io.File)
	 */
	@Override
	public DeBruijnSubgraphGexfExporter saveTo(File file) {
		GexfHelper.saveTo(gexf, file);
		return this;
	}
}
