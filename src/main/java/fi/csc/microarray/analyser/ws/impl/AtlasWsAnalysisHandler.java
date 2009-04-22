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
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import fi.csc.microarray.analyser.ws.HtmlUtil;
import fi.csc.microarray.analyser.ws.ResultTableCollector;
import fi.csc.microarray.util.XmlUtil;

public class AtlasWsAnalysisHandler {


	public static void main(String[] args) throws SAXException, ParserConfigurationException, TransformerException, SOAPException, IOException {
		execute(new String[] {"ENSG00000125304"});
	}

	private static void execute(String[] probes) throws SAXException, ParserConfigurationException, TransformerException, SOAPException, IOException {
		ResultTableCollector annotations = query(probes);
		HtmlUtil.writeHtmlTable(annotations, new String[] {"updn", "experiment_accession", "experiment_description", "gene_name"});
	}
	
	public static ResultTableCollector query(String[] genes) throws SAXException, ParserConfigurationException, TransformerException, SOAPException, IOException {
		try {
			ResultTableCollector annotationCollector = new ResultTableCollector();
			MessageFactory mf = MessageFactory.newInstance();
			SOAPMessage request = mf.createMessage();
			SOAPPart soapPart = request.getSOAPPart();
			SOAPEnvelope soapEnvelope = soapPart.getEnvelope();
			soapEnvelope.addNamespaceDeclaration("web", "http://webservices.service.ae3/");
			SOAPBody soapBody = soapEnvelope.getBody();
			
			SOAPElement batchQueryElement = soapBody.addChildElement("batchQuery", "web");
			
			SOAPElement genesElement = batchQueryElement.addChildElement("q_genes", "web");
			
			for (String gene : genes) { 
				SOAPElement string = genesElement.addChildElement("string", "web");
				string.setTextContent(gene);
			}

			batchQueryElement.addChildElement("q_expts", "web");
			batchQueryElement.addChildElement("q_orgn", "web");
			batchQueryElement.addChildElement("q_updn", "web");

			request.saveChanges();
			SOAPConnectionFactory connectionFactory = SOAPConnectionFactory.newInstance();
			SOAPConnection soapConnection = connectionFactory.createConnection();
			URL endpoint = new URL("http://www.ebi.ac.uk/microarray-as/atlas/services/AtlasWebService");

			SOAPMessage resp = soapConnection.call(request, endpoint);
			
			ByteArrayOutputStream out = new ByteArrayOutputStream();
			resp.writeTo(out);
			soapConnection.close();
			
			Document response = XmlUtil.getInstance().parseReader(new InputStreamReader(new ByteArrayInputStream(out.toByteArray())));

			NodeList rows = response.getDocumentElement().getElementsByTagName("ns1:anyType2anyTypeMap");
			
			for (int i = 0; i < rows.getLength(); i++) {
				Element row = (Element)rows.item(i);
				NodeList fields = row.getElementsByTagName("ns1:entry");
				for (int j = 0; j < fields.getLength(); j++) {
					Element field = (Element)fields.item(j);
					String name = field.getChildNodes().item(0).getTextContent();
					String value = field.getChildNodes().item(1).getTextContent();
					annotationCollector.addAnnotation(i, name, value);
				}
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
