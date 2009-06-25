package fi.csc.microarray.analyser.ws.impl;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.soap.MessageFactory;
import javax.xml.soap.SOAPBody;
import javax.xml.soap.SOAPConnection;
import javax.xml.soap.SOAPConnectionFactory;
import javax.xml.soap.SOAPElement;
import javax.xml.soap.SOAPEnvelope;
import javax.xml.soap.SOAPException;
import javax.xml.soap.SOAPMessage;
import javax.xml.soap.SOAPPart;
import javax.xml.transform.TransformerException;

import org.w3c.dom.Document;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import fi.csc.microarray.analyser.ws.HtmlUtil;
import fi.csc.microarray.analyser.ws.ResultTableCollector;
import fi.csc.microarray.analyser.ws.ResultTableCollector.ResultRow;
import fi.csc.microarray.analyser.ws.ResultTableCollector.RowFilter;
import fi.csc.microarray.util.XmlUtil;

public class CpdbWsUtils {


	public static void main(String[] args) throws SAXException, ParserConfigurationException, TransformerException, SOAPException, IOException {
		String[] probes = new String[] {"TM9SF2", "FOLR3", "IER2", "TMED2", "TMEM131", "PVRL2", "MIA3"};
		ResultTableCollector annotations = query(probes);
		annotations.filterRows(new RowFilter() {
			public boolean shouldRemove(ResultRow row) {
				return Double.parseDouble(row.getValue("ns1:pValue")) > 0.0005d;
			}
		});
		HtmlUtil.writeHtmlTable(annotations, new String[] {"ns1:pValue", "ns1:pathway", "ns1:database"}, "ConsensusPathDB annotation", System.out);
	}

	public static ResultTableCollector query(String[] genes) throws SAXException, ParserConfigurationException, TransformerException, SOAPException, IOException {
		try {
			ResultTableCollector annotationCollector = new ResultTableCollector();
			MessageFactory mf = MessageFactory.newInstance();
			SOAPMessage request = mf.createMessage();
			SOAPPart soapPart = request.getSOAPPart();
			SOAPEnvelope soapEnvelope = soapPart.getEnvelope();
			soapEnvelope.addNamespaceDeclaration("cpdb", "http://cpdb.molgen.mpg.de/cpdb_1_04");
			SOAPBody soapBody = soapEnvelope.getBody();
			
			SOAPElement elementORApathways = soapBody.addChildElement("ORApathways", "cpdb");
			
			SOAPElement elementIdType = elementORApathways.addChildElement("CPDB_idType", "cpdb");
			elementIdType.setTextContent("hgnc");
			
			for (String gene : genes) { 
				SOAPElement elementIdList1 = elementORApathways.addChildElement("CPDB_idList", "cpdb");
				elementIdList1.setTextContent(gene);
			}

			SOAPElement elementBol = elementORApathways.addChildElement("bol", "cpdb");
			elementBol.setTextContent("1");

			request.saveChanges();
			SOAPConnectionFactory connectionFactory = SOAPConnectionFactory.newInstance();
			SOAPConnection soapConnection = connectionFactory.createConnection();
			URL endpoint = new URL("http://cpdb.molgen.mpg.de/soap");

			SOAPMessage resp = soapConnection.call(request, endpoint);
			
			ByteArrayOutputStream out = new ByteArrayOutputStream();
			resp.writeTo(out);
			soapConnection.close();
			
			Document response = XmlUtil.getInstance().parseReader(new InputStreamReader(new ByteArrayInputStream(out.toByteArray())));
			
			NodeList childNodes = response.getDocumentElement().getChildNodes().item(1).getChildNodes().item(0).getChildNodes();;

			int index = 0;
			String prevName = null; 
			for (int i = 0; i < childNodes.getLength(); i++) {
				String name = childNodes.item(i).getNodeName();
				String value = childNodes.item(i).getTextContent();

				if (name.equals(prevName)) {
					index++;
				} else {
					index = 0; // new field
				}

				annotationCollector.addAnnotation(index, name, value);
				prevName = name;
			}
			
			return annotationCollector;

		} catch (java.io.IOException ioe) {
			ioe.printStackTrace();
			throw ioe;
		} catch (SOAPException soape) {
			soape.printStackTrace();
			throw soape;
		}
	}


}
