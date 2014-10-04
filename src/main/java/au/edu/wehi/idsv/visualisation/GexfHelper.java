package au.edu.wehi.idsv.visualisation;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;
import it.uniroma1.dis.wsngroup.gexf4j.core.Gexf;
import it.uniroma1.dis.wsngroup.gexf4j.core.GexfWriter;
import it.uniroma1.dis.wsngroup.gexf4j.core.IntervalType;
import it.uniroma1.dis.wsngroup.gexf4j.core.Node;
import it.uniroma1.dis.wsngroup.gexf4j.core.data.Attribute;
import it.uniroma1.dis.wsngroup.gexf4j.core.data.AttributeValue;
import it.uniroma1.dis.wsngroup.gexf4j.core.data.AttributeValueList;
import it.uniroma1.dis.wsngroup.gexf4j.core.impl.StaxGraphWriter;

public class GexfHelper {
	private static final Log log = Log.getInstance(GexfHelper.class);
	public static void setDynamicAttribute(Node node, int timeStamp, Attribute attr, Object value) {
		AttributeValueList attrValues = node.getAttributeValues();
		AttributeValue last = getLatestDynamicAttribute(attrValues, attr);
		assert(last == null || (int)(Integer)last.getStartValue() <= timeStamp);
		if (last == null) {
			if (value != null) {
				// create new attribute
				attrValues.createValue(attr, value.toString())
				.setStartIntervalType(IntervalType.CLOSE)
				.setStartValue(timeStamp);
			}
		} else {
			if (value == null) {
				// attribute no longer applies
				last.setEndValue(timeStamp);
				last.setEndIntervalType(IntervalType.OPEN);
			} else {
				String newValue = value.toString();
				if (last.getValue().equals(newValue)) {
					// attribute value has not changed
				} else if ((int)(Integer)last.getStartValue() == timeStamp) {
					// overwrite with our new value
					last.setValue(newValue);
				} else {
					// new attribute value
					last.setEndValue(timeStamp);
					last.setEndIntervalType(IntervalType.OPEN);
					last = attrValues.createValue(attr, newValue)
						.setStartIntervalType(IntervalType.CLOSE)
						.setStartValue(timeStamp);
				}
				last.clearEndDate();
			}
		}
	}
	public static AttributeValue getLatestDynamicAttribute(AttributeValueList attrValues, Attribute attr) {
		AttributeValue best = null;
		for (int i = 0 ; i < attrValues.size(); i++) {
			AttributeValue av = attrValues.get(i);
			if (av.getAttribute() == attr) {
				if (best == null || (int)(Integer)best.getStartValue() < (int)(Integer)av.getStartValue()) {
					best = av;
				}
			}
		}
		return best;
	}
	public static void saveTo(Gexf graph, File file) {
		file.getParentFile().mkdir();
		GexfWriter graphWriter = new StaxGraphWriter();
		Writer out = null;
		try {
			out = new FileWriter(file, false);
			graphWriter.writeToStream(graph, out, "UTF-8");
		} catch (IOException e) {
			log.error(String.format("Error writing graph visualisation to %s", file), e);
		} finally {
			CloserUtil.close(out);
			CloserUtil.close(graphWriter);
		}
	}
}
