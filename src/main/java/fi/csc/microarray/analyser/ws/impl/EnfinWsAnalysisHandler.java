package fi.csc.microarray.analyser.ws.impl;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.HashMap;

import javax.xml.namespace.QName;
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

public class EnfinWsAnalysisHandler {


	public static void main(String[] args) throws SAXException, ParserConfigurationException, TransformerException, SOAPException, IOException {
		execute(new String[] {"P38398", "P38398"});
	}

	private static void execute(String[] probes) throws SAXException, ParserConfigurationException, TransformerException, SOAPException, IOException {
		ResultTableCollector annotations = query(probes);
		HtmlUtil.writeHtmlTable(annotations, new String[] {"Interaction", "Participants"}, "Enfin IntAct annotation", new File("test.html"));
	}

	public static ResultTableCollector query(String[] proteins) throws SAXException, ParserConfigurationException, TransformerException, SOAPException, IOException {
		try {
			ResultTableCollector annotationCollector = new ResultTableCollector();
			MessageFactory mf = MessageFactory.newInstance();
			SOAPMessage request = mf.createMessage();
			SOAPPart soapPart = request.getSOAPPart();
			SOAPEnvelope soapEnvelope = soapPart.getEnvelope();
			soapEnvelope.addNamespaceDeclaration("int", "http://ebi.ac.uk/enfin/core/web/services/intact");
			soapEnvelope.addNamespaceDeclaration("model", "http://ebi.ac.uk/enfin/core/model");
			
			SOAPBody soapBody = soapEnvelope.getBody();
			
			SOAPElement elementFindPartners = soapBody.addChildElement("findPartners", "int");
			
			SOAPElement elementEntries = elementFindPartners.addChildElement("entries", "model");
			
			SOAPElement elementEntry = elementEntries.addChildElement("entry", "model");
			
			int id = 1;
			for (String protein : proteins) { 
				SOAPElement elementMolecule = elementEntry.addChildElement("molecule", "model");
				elementMolecule.addAttribute(new QName("id"), "ID" + id);
				
				SOAPElement elementXrefs = elementMolecule.addChildElement("xrefs", "model");
				SOAPElement elementPrimaryRef = elementXrefs.addChildElement("primaryRef", "model");
				elementPrimaryRef.addAttribute(new QName("refTypeAc"), "MI:0358");
				elementPrimaryRef.addAttribute(new QName("refType"), "primary-reference");
				elementPrimaryRef.addAttribute(new QName("id"), protein);
				elementPrimaryRef.addAttribute(new QName("db"), "not specified");

				SOAPElement elementMoleculeType = elementMolecule.addChildElement("moleculeType", "model");
				elementMoleculeType.addAttribute(new QName("termAc"), "MI:0326");
				elementMoleculeType.addAttribute(new QName("term"), "protein");

				id++;
			}

			SOAPElement elementSet = elementEntry.addChildElement("set", "model");
			elementSet.addAttribute(new QName("id"), "ID" + id);
			for (int i = 1; i < id; i++) {
				SOAPElement elementParticipant = elementSet.addChildElement("participant", "model");
				elementParticipant.addAttribute(new QName("moleculeRef"), "ID" + i);
			}

			id++;
			SOAPElement elementExperiment = elementEntry.addChildElement("experiment", "model");
			elementEntry.addAttribute(new QName("id"), "ID" + id);

			SOAPElement elementNames = elementExperiment.addChildElement("names", "model");
			SOAPElement elementShortLabel = elementNames.addChildElement("shortLabel", "model");
			elementShortLabel.setTextContent("enfin-utility");
			SOAPElement elementFullName = elementNames.addChildElement("fullName", "model");
			elementFullName.setTextContent("Utility service to create Enfin XML document from given simple input (expected: various protein IDs).");
			SOAPElement elementResult = elementExperiment.addChildElement("result", "model");
			elementResult.setTextContent("ID" + (id - 1));
         
			id++;
			SOAPElement elementParameter = elementEntry.addChildElement("parameter", "model");
			elementParameter.addAttribute(new QName("factor"), "" + id);
			elementParameter.addAttribute(new QName("term"), "IdCounter");

			request.saveChanges();
			SOAPConnectionFactory connectionFactory = SOAPConnectionFactory.newInstance();
			SOAPConnection soapConnection = connectionFactory.createConnection();
			URL endpoint = new URL("http://www.ebi.ac.uk/enfin-srv/encore/intact/service");

			SOAPMessage resp = soapConnection.call(request, endpoint);
			
			ByteArrayOutputStream out = new ByteArrayOutputStream();
			resp.writeTo(out);
			soapConnection.close();
			
			Document response = XmlUtil.getInstance().parseReader(new InputStreamReader(new ByteArrayInputStream(out.toByteArray())));

			NodeList childNodes = response.getDocumentElement().getChildNodes().item(0).getChildNodes().item(0).getChildNodes().item(0).getChildNodes().item(0).getChildNodes();

			// iterate over molecules
			HashMap<String, String> moleculeMap = new HashMap<String, String>();
			for (int i = 0; i < childNodes.getLength(); i++) {
				if ("molecule".equals(childNodes.item(i).getNodeName())) {
					Element molecule = (Element)childNodes.item(i);					
					String moleculeId = molecule.getAttribute("id");
					
					Element primaryRef = (Element)molecule.getElementsByTagName("xrefs").item(0).getChildNodes().item(0);
					String moleculeName = primaryRef.getAttribute("id");
					moleculeMap.put(moleculeId, moleculeName);		
				}
			}

			// iterate over interactions
			int index = 0;
			for (int i = 0; i < childNodes.getLength(); i++) {
				if ("set".equals(childNodes.item(i).getNodeName())) {
					Element set = (Element)childNodes.item(i);
					NodeList names = set.getElementsByTagName("names");
					if (names.getLength() > 0 && "IntAct interaction".equals(names.item(0).getChildNodes().item(0).getTextContent())) {
						Element primaryRef = (Element)set.getElementsByTagName("xrefs").item(0).getChildNodes().item(0);						
						annotationCollector.addAnnotation(index, "Interaction", primaryRef.getAttribute("id"));
						
						String participantValue = "";
						NodeList participants = set.getElementsByTagName("participant");
						for (int p = 0; p < participants.getLength(); p++) {
							Element participant = (Element)participants.item(p);
							String moleculeName = moleculeMap.get(participant.getAttribute("moleculeRef"));
							participantValue += (" " + moleculeName);
						}
						annotationCollector.addAnnotation(index, "Participants", participantValue);
						
						index++;
					}						
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
