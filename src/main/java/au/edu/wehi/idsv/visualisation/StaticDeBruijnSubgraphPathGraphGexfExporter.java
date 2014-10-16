package au.edu.wehi.idsv.visualisation;

import it.uniroma1.dis.wsngroup.gexf4j.core.Node;
import it.uniroma1.dis.wsngroup.gexf4j.core.data.Attribute;
import it.uniroma1.dis.wsngroup.gexf4j.core.data.AttributeList;
import it.uniroma1.dis.wsngroup.gexf4j.core.data.AttributeType;
import au.edu.wehi.idsv.debruijn.subgraph.DeBruijnSubgraphNode;
import au.edu.wehi.idsv.debruijn.subgraph.SubgraphPathNode;

public class StaticDeBruijnSubgraphPathGraphGexfExporter extends StaticDeBruijnPathGraphGexfExporter<DeBruijnSubgraphNode, SubgraphPathNode> {
	private final Attribute attrIsReference;
	private final Attribute attrIsNonReference;
	public StaticDeBruijnSubgraphPathGraphGexfExporter(int k) {
		super(k);
		AttributeList nodeStaticAttrList = getNodeStaticAttributeList();
		attrIsReference = nodeStaticAttrList.createAttribute("r", AttributeType.BOOLEAN, "containsReferenceKmer");
		attrIsNonReference = nodeStaticAttrList.createAttribute("nr", AttributeType.BOOLEAN, "containsNonReferenceKmer");
	}
	@Override
	protected void setStaticAttributes(Node node, SubgraphPathNode pn) {
		super.setStaticAttributes(node, pn);
		node.getAttributeValues().createValue(attrIsReference, ((Boolean)pn.containsReferenceKmer()).toString());
		node.getAttributeValues().createValue(attrIsNonReference, ((Boolean)pn.containsNonReferenceKmer()).toString());
	}
}
