package fi.csc.microarray.analyser.ws.impl;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

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
import fi.csc.microarray.analyser.ws.HtmlUtil.ValueHtmlFormatter;
import fi.csc.microarray.util.Strings;
import fi.csc.microarray.util.XmlUtil;

public class EnfinWsUtils {

	static interface AnnotationIdentifier {
		boolean isAnnotation(Element setElement);
	}

	static interface AnnotationNameFinder {
		String findAnnotationName(Element setElement);
	}

	public static void main(String[] args) throws SAXException, ParserConfigurationException, TransformerException, SOAPException, IOException {
		String[] probes = new String[] {		
				"204704_s_at",
				"221589_s_at",
				"206065_s_at",
				"209459_s_at",
				"209460_at",
				"206024_at",
				"205719_s_at",
				"205892_s_at",
				"202036_s_at",
				"206054_at",
				"209443_at"
		};
//		String[] probes = JavaJobUtils.getProbes(new File("/tmp/two-sample.tsv"));
		
		ResultTableCollector intactAnnotations = queryIntact(probes);
		writeIntactResult(intactAnnotations, new File("intact.html"), new File("intact.tsv"));
		
		ResultTableCollector reactomeAnnotations = queryReactome(probes);
		writeReactomeResult(reactomeAnnotations, new File("reactome.html"), new File("reactome.tsv"));
	}

	public static void writeIntactResult(ResultTableCollector intactAnnotations, File htmlFile, File textFile) throws FileNotFoundException {
		ValueHtmlFormatter interactionIdFormatter = new ValueHtmlFormatter() {
			public String format(String string, String[] currentRow) {
				return "<a href=\"http://www.ebi.ac.uk/intact/pages/interactions/interactions.xhtml?conversationContext=1&queryTxt=" + string.replace(' ', '+') + "\">" + string + "</a>";
			}
		};
		
		ValueHtmlFormatter uniprotFormatter = new ValueHtmlFormatter() {
			public String format(String string, String[] currentRow) {
				String formatted = "";
				for (String protein : string.split(" ")) {
					formatted += "<a href=\"http://www.uniprot.org/uniprot/" + protein.replace(' ', '+') + "\">" + protein + "</a> "; 
				}
				return formatted;
			}
		};
		
		HtmlUtil.writeHtmlTable(intactAnnotations, new String[] {"Name", "Probe IDs", "Participants"}, new String[] {"Interaction", "Probe ID", "Interacting proteins"},  new HtmlUtil.ValueHtmlFormatter[] {interactionIdFormatter, HtmlUtil.NO_FORMATTING_FORMATTER, uniprotFormatter}, "IntAct protein interactions", new FileOutputStream(htmlFile));
		HtmlUtil.writeTextTable(intactAnnotations, new String[] {"Name", "Probe IDs", "Participants"}, new String[] {"Interaction", "Probe ID", "Interacting proteins"}, new FileOutputStream(textFile));
	}

	public static void writeReactomeResult(ResultTableCollector reactomeAnnotations, File htmlFile, File textFile) throws FileNotFoundException {
		ValueHtmlFormatter pathwayNameFormatter = new ValueHtmlFormatter() {
			public String format(String string, String[] currentRow) {
				return "<a href=\"http://www.reactome.org/cgi-bin/search2?DB=gk_current&OPERATOR=ALL&QUERY=" + string.replace(' ', '+') + "&SPECIES=&SUBMIT=Go!\">" + string + "</a>";
			}
		};
		
		HtmlUtil.writeHtmlTable(reactomeAnnotations, new String[] { "Name", "Probe IDs" }, new String[] {"Pathway", "Probe ID's for participating proteins"}, new HtmlUtil.ValueHtmlFormatter[] {pathwayNameFormatter, HtmlUtil.NO_FORMATTING_FORMATTER}, "Reactome pathway associations", new FileOutputStream(htmlFile));
		HtmlUtil.writeTextTable(reactomeAnnotations, new String[] { "Name", "Probe IDs" }, new String[] {"Pathway", "Probe ID's for participating proteins"}, new FileOutputStream(textFile));
	}

	public static ResultTableCollector queryIntact(String[] probes) throws SOAPException, MalformedURLException, SAXException, IOException, ParserConfigurationException, TransformerException {
		
		Document uniprotResponse = queryUniprotIds(probes);

		// query IntAct with UniProt identifiers
		SOAPMessage intactSoapMessage = initialiseSoapMessage();
		SOAPBody intactSoapBody = initialiseSoapBody(intactSoapMessage);
		
		attachEnfinXml(intactSoapBody, fetchEnfinXml(uniprotResponse), "findPartners", "http://ebi.ac.uk/enfin/core/web/services/intact");
		
		Document intactResponse = sendSoapMessage(intactSoapMessage, new URL("http://www.ebi.ac.uk/enfin-srv/encore/intact/service"));
		
		return collectAnnotations(intactResponse, new AnnotationIdentifier() {

			public boolean isAnnotation(Element setElement) {
				List<Element> names = XmlUtil.getChildElements(setElement, "names");
				return !names.isEmpty() && "IntAct interaction".equals(XmlUtil.getChildElement(names.get(0), "fullName").getTextContent());
			}
			
		}, new AnnotationNameFinder() {

			public String findAnnotationName(Element setElement) {
				Element primaryRef = (Element)XmlUtil.getChildElement(setElement, "xrefs").getChildNodes().item(0);
				return primaryRef.getAttribute("id");
			}
			
		});
	}

	public static ResultTableCollector queryReactome(String[] probes) throws SOAPException, MalformedURLException, SAXException, IOException, ParserConfigurationException, TransformerException {
		
		Document uniprotResponse = queryUniprotIds(probes);

		// query Reactome with UniProt identifiers
		SOAPMessage intactSoapMessage = initialiseSoapMessage();
		SOAPBody intactSoapBody = initialiseSoapBody(intactSoapMessage);
		
		attachEnfinXml(intactSoapBody, fetchEnfinXml(uniprotResponse), "findPath", "http://ebi.ac.uk/enfin/core/web/services/reactome");
		
		Document reactomeResponse = sendSoapMessage(intactSoapMessage, new URL("http://www.ebi.ac.uk/enfin-srv/encore/reactome/service"));
//		XmlUtil.printXml(reactomeResponse, System.out);
		
		return collectAnnotations(reactomeResponse, new AnnotationIdentifier() {

			public boolean isAnnotation(Element setElement) {
				List<Element> setTypes = XmlUtil.getChildElements(setElement, "setType");
				return !setTypes.isEmpty() && "Reactome".equals(setTypes.get(0).getAttribute("db"));
			}
			
		}, new AnnotationNameFinder() {

			public String findAnnotationName(Element setElement) {
				Element fullName = (Element)XmlUtil.getChildElement(setElement, "names").getChildNodes().item(0);
				return fullName.getTextContent();
			}
			
		});
	}

	private static Document queryUniprotIds(String[] probes) throws SOAPException, SAXException, IOException, ParserConfigurationException, MalformedURLException {
		
		// Step 1. Create ENFIN XML out of Affy probe list
		SOAPMessage probeSoapMessage = initialiseSoapMessage();
		SOAPBody probeSoapBody = initialiseSoapBody(probeSoapMessage);

		SOAPElement probeEntries = createOperation(probeSoapBody, "docFromAffyList", "http://ebi.ac.uk/enfin/core/web/services/utility");

		for (String probe : probes) {
			SOAPElement arg = probeEntries.addChildElement("parameter");
			arg.setTextContent(probe);
		}
		
		Document probeResponse = sendSoapMessage(probeSoapMessage, new URL("http://www.ebi.ac.uk/enfin-srv/encore/utility/service"));
		
		// Step 2. Convert Affy probes to UniProt identifiers
		SOAPMessage uniprotSoapMessage = initialiseSoapMessage();
		SOAPBody uniprotSoapBody = initialiseSoapBody(uniprotSoapMessage);
		
		attachEnfinXml(uniprotSoapBody, fetchEnfinXml(probeResponse), "mapAffy2UniProt", "http://ebi.ac.uk/enfin/core/web/services/affy2uniprot");

		Document uniprotResponse = sendSoapMessage(uniprotSoapMessage, new URL("http://www.ebi.ac.uk/enfin-srv/encore/affy2uniprot/service"));
		return uniprotResponse;
	}

	private static Document fetchEnfinXml(Document response) throws ParserConfigurationException {
		Document document = XmlUtil.newDocument();
		document.appendChild(document.importNode(response.getDocumentElement().getElementsByTagNameNS("http://ebi.ac.uk/enfin/core/model", "entries").item(0), true));
		return document;
	}
	
	private static ResultTableCollector collectAnnotations(Document response, AnnotationIdentifier annotationIdentifier, AnnotationNameFinder annotationNameFinder) {
		ResultTableCollector annotationCollector = new ResultTableCollector();
		NodeList childNodes = response.getDocumentElement().getChildNodes().item(0).getChildNodes().item(0).getChildNodes().item(0).getChildNodes().item(0).getChildNodes();

		// iterate over molecules
		HashMap<String, String> moleculeMap = new HashMap<String, String>();
		for (int i = 0; i < childNodes.getLength(); i++) {
			if ("molecule".equals(childNodes.item(i).getLocalName())) {
				Element molecule = (Element)childNodes.item(i);					
				String moleculeId = molecule.getAttribute("id");
				
				Element primaryRef = (Element)XmlUtil.getChildElement(molecule, "xrefs").getChildNodes().item(0);
				String moleculeName = primaryRef.getAttribute("id");
				moleculeMap.put(moleculeId, moleculeName);		
			}
		}

		// iterate over mappings (Affy->UniProt)
		HashMap<String, String> proteinToAffyMap = new HashMap<String, String>();
		for (int i = 0; i < childNodes.getLength(); i++) {
			if ("set".equals(childNodes.item(i).getLocalName())) {
				Element set = (Element) childNodes.item(i);
				Element setType = XmlUtil.getChildElement(set, "setType");
				if (setType != null && "Affymetrix ID mapped to UniProt accession(s)".equals(setType.getAttribute("term"))) {
					List<Element> participants = XmlUtil.getChildElements(set, "participant");
					String affyRef = moleculeMap.get(participants.get(0).getAttribute("moleculeRef"));
					// iterate over rest
					for (int j = 1; j < participants.size(); j++) {
						String proteinId = participants.get(j).getAttribute("moleculeRef");
						proteinToAffyMap.put(proteinId, affyRef);
					}
				}
			}
		}
		
		// iterate over interactions
		int index = 0;
		for (int i = 0; i < childNodes.getLength(); i++) {
			if ("set".equals(((Element)childNodes.item(i)).getLocalName())) {
				Element set = (Element)childNodes.item(i);
				
				if (annotationIdentifier.isAnnotation(set)) {
					
					String annotationName = annotationNameFinder.findAnnotationName(set);
					annotationCollector.addAnnotation(index, "Name", annotationName);
					
					Set<String> probeids = new HashSet<String>();
					Set<String> moleculeNames = new HashSet<String>();
					List<Element> participants = XmlUtil.getChildElements(set, "participant");
					for (Element participant : participants) {
						String participantValue = participant.getAttribute("moleculeRef");
						String moleculeName = moleculeMap.get(participantValue);
						moleculeNames.add(moleculeName);
						String probeName = proteinToAffyMap.get(participantValue);
						if (probeName != null) {
							probeids.add(probeName);
						}
					}
					annotationCollector.addAnnotation(index, "Probe IDs", Strings.delimit(probeids, " "));
					annotationCollector.addAnnotation(index, "Participants", Strings.delimit(moleculeNames, " "));
					
					index++;
				}						
			}
		}
		
		return annotationCollector;

	}
	

	private static Document sendSoapMessage(SOAPMessage soapMessage, URL endpoint) throws SAXException, IOException, ParserConfigurationException, SOAPException {
		soapMessage.saveChanges();
		SOAPConnectionFactory connectionFactory = SOAPConnectionFactory.newInstance();
		SOAPConnection soapConnection = connectionFactory.createConnection();

		SOAPMessage resp = soapConnection.call(soapMessage, endpoint);

		ByteArrayOutputStream out = new ByteArrayOutputStream();
		resp.writeTo(out);
		soapConnection.close();
		
		return XmlUtil.parseReader(new InputStreamReader(new ByteArrayInputStream(out.toByteArray())));
	}

	private static SOAPBody initialiseSoapBody(SOAPMessage message) throws SOAPException {
		SOAPPart soapPart = message.getSOAPPart();
		SOAPEnvelope soapEnvelope = soapPart.getEnvelope();
		
		
		return soapEnvelope.getBody();
	}
	
	private static SOAPMessage initialiseSoapMessage() throws SOAPException {
		MessageFactory mf = MessageFactory.newInstance();
		return mf.createMessage();
	}
	
	private static void attachEnfinXml(SOAPBody soapBody, Document enfinXml, String operation, String operationNamespace) throws SOAPException, ParserConfigurationException {
		Document document = XmlUtil.newDocument();
		document.appendChild(document.createElementNS(operationNamespace, operation));	
		document.getDocumentElement().appendChild(document.importNode(enfinXml.getDocumentElement(), true));		
		soapBody.addDocument(document);
	}

	private static SOAPElement createOperation(SOAPBody soapBody, String operation, String operationNamespace) throws SOAPException {
		soapBody.addNamespaceDeclaration("oper", operationNamespace);
		SOAPElement elementFindPartners = soapBody.addChildElement(operation, "oper");
		elementFindPartners.addNamespaceDeclaration("model", "http://ebi.ac.uk/enfin/core/model");
		return elementFindPartners;
	}
}
