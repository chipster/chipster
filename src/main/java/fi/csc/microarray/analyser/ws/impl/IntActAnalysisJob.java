package fi.csc.microarray.analyser.ws.impl;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.HashMap;

import javax.xml.soap.MessageFactory;
import javax.xml.soap.SOAPBody;
import javax.xml.soap.SOAPConnection;
import javax.xml.soap.SOAPConnectionFactory;
import javax.xml.soap.SOAPElement;
import javax.xml.soap.SOAPEnvelope;
import javax.xml.soap.SOAPMessage;
import javax.xml.soap.SOAPPart;

import org.w3c.dom.Document;
import org.w3c.dom.NodeList;

import fi.csc.microarray.analyser.ws.SoapAnalysisJob;
import fi.csc.microarray.analyser.ws.resultgraph.ResultGraph;
import fi.csc.microarray.util.XmlUtil;

public class IntActAnalysisJob extends SoapAnalysisJob {

	@Override
	public Document getEnvelope(Iterable<String> query) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Document getResultTable(ResultGraph resultGraph) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String getVVSADL() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	protected void execute() throws Exception {
		MessageFactory mf = MessageFactory.newInstance();
		SOAPMessage soapMessage = mf.createMessage();
		SOAPPart soapPart = soapMessage.getSOAPPart();
		SOAPEnvelope soapEnvelope = soapPart.getEnvelope();
		soapEnvelope.addNamespaceDeclaration("cpdb", "http://cpdb.molgen.mpg.de/cpdb_1_04");
		SOAPBody soapBody = soapEnvelope.getBody();
		
		SOAPElement elementORApathways = soapBody.addChildElement("ORApathways", "cpdb");
		
		SOAPElement elementIdType = elementORApathways.addChildElement("CPDB_idType", "cpdb");
		elementIdType.setTextContent("hgnc");
		
		SOAPElement elementIdList= elementORApathways.addChildElement("CPDB_idList", "cpdb");
		elementIdList.setTextContent("annamulletuloksia");

		SOAPElement elementBol = elementORApathways.addChildElement("CPDB_bol", "cpdb");
		elementBol.setTextContent("1");

		soapMessage.saveChanges();
		SOAPConnectionFactory connectionFactory = SOAPConnectionFactory.newInstance();
		SOAPConnection soapConnection = connectionFactory.createConnection();
		URL endpoint = new URL("http://cpdb.molgen.mpg.de/soap");

		SOAPMessage resp = soapConnection.call(soapMessage, endpoint);
		
		ByteArrayOutputStream out = new ByteArrayOutputStream();
		resp.writeTo(out);
		soapConnection.close();
		
		Document response = XmlUtil.getInstance().parseReader(new InputStreamReader(new ByteArrayInputStream(out.toByteArray())));
		
		NodeList childNodes = response.getDocumentElement().getChildNodes().item(1).getChildNodes().item(0).getChildNodes();;
		HashMap<String, Integer> counter = new HashMap<String, Integer>();
		counter.put("ns1:pValue", 0);
		counter.put("ns1:pathway", 0);
		counter.put("ns1:database", 0);
		counter.put("ns1:pathway_members", 0);
		counter.put("ns1:membersInputOverlap", 0);			
		counter.put("ns1:pathwaySize", 0);
		counter.put("ns1:overlapSize", 0);
		
		for (int i = 0; i < childNodes.getLength(); i++) {
			String name = childNodes.item(i).getNodeName();
			counter.put(name, counter.get(name).intValue() + 1);
		}
		
		System.out.println(counter);		
		
	}

}
