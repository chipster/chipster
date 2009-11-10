package fi.csc.microarray.util;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.Reader;
import java.io.UnsupportedEncodingException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.List;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

/**
 *
 * @author  Aleksi Kallio
 */
public class XmlUtil {
    
	/**
	 * Should be deprecated some day.
	 */
    public static synchronized XmlUtil getInstance() {
        return new XmlUtil();
    }
        
    public static Document newDocument() throws ParserConfigurationException {
        return newDocumentBuilder().newDocument();
    }
    
    public static Document parseFile(File file) throws org.xml.sax.SAXException, IOException, ParserConfigurationException {
        return newDocumentBuilder().parse(file);
    }
    
    public static Document parseReader(Reader reader) throws org.xml.sax.SAXException, IOException, ParserConfigurationException {
        return newDocumentBuilder().parse(new org.xml.sax.InputSource(reader));
    }
    
    public static void printXml(Document xml, Writer out) throws TransformerException, UnsupportedEncodingException {
        Transformer transformer = TransformerFactory.newInstance().newTransformer();
        transformer.setOutputProperty(OutputKeys.METHOD, "xml");
        transformer.setOutputProperty(OutputKeys.INDENT, "yes");        
        
        transformer.transform(new DOMSource(xml), new StreamResult(out));        
    }
    
	public static Element getChildWithAttributeValue(Element parent, String attrName, String attrValue) {
		NodeList childNodes = parent.getChildNodes();
		for (int i = 0; i < childNodes.getLength(); i++) {
			Node node = childNodes.item(i);
			if (node instanceof Element) {
				Element element = (Element)node;
				if (attrValue.equals(element.getAttribute(attrName))) {
					return element;
				}
			}
		}
		return null;
	}

	public static Element getChildWithAttribute(Element parent, String attrName) {
		NodeList childNodes = parent.getChildNodes();
		for (int i = 0; i < childNodes.getLength(); i++) {
			Node node = childNodes.item(i);
			if (node instanceof Element) {
				Element element = (Element)node;
				if (!"".equals(element.getAttribute(attrName))) {
					return element;
				}
			}
		}
		return null;
	}

	public static void printXml(Document response, OutputStream out) throws UnsupportedEncodingException, TransformerException {
		printXml(response, new OutputStreamWriter(out));		
	}

	/**
	 * Convenience method for getting all child elements.
	 * 
	 * @see #getChildElements(Element, String)
	 */
	public static List<Element> getChildElements(Element parent) {
		return getChildElements(parent, null);
	}
	
	/**
	 * Gets the child elements of a parent element. Unlike DOM's getElementsByTagName, this does no recursion,
	 * uses local name (namespace free) instead of tag name, result is a proper Java data structure and result
	 * needs no casting. In other words, this method does not suck unlike DOM.
	 * 
	 * @param parent the XML parent element
	 * @param name name of the child elements, if null then all are returned
	 */
	public static List<Element> getChildElements(Element parent, String name) {
		List<Element> childElements = new ArrayList<Element>();
		NodeList childNodes = parent.getChildNodes();
		
		for (int i = 0; i < childNodes.getLength(); i++) {
			// get elements
			if (childNodes.item(i).getNodeType() == Node.ELEMENT_NODE) {
				
				// match element name
				Element childElement = (Element) childNodes.item(i);
				if (name == null || childElement.getLocalName().equals(name)) {
					childElements.add(childElement);
				}
			}
		}
		
		return childElements;
	}

	public static Element getChildElement(Element parent, String name) {
		return getChildElement(parent, name, false);
	}
	
	public static Element getChildElement(Element parent, String name, boolean strict) {
		List<Element> childElements = getChildElements(parent, name);
		if (strict && childElements.size() != 1) {
			throw new IllegalArgumentException("parent must contain exactly one element with the given name");
		} 
		
		return childElements.isEmpty() ? null : childElements.get(0);	
	}
	
    private static DocumentBuilder newDocumentBuilder() throws ParserConfigurationException {
    	// SAXParsers are not concurrency compatible, so always return a new instance to prevent thread issues 
        DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
        dbf.setNamespaceAware(true);
        return dbf.newDocumentBuilder();
    }
    

}
