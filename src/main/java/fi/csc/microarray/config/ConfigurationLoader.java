package fi.csc.microarray.config;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import fi.csc.microarray.util.XmlUtil;

public class ConfigurationLoader {
	
	private static final String VERSION_ATTRIBUTE = "content-version";
	
	public static class OldConfigurationFormatException extends Exception {

		public OldConfigurationFormatException(String string) {
			super(string);
		}
		
	}
	public static void addFromFile(ConfigurationModule configuration, File file, int requiredVersion) throws IOException, SAXException, ParserConfigurationException, OldConfigurationFormatException {
		FileInputStream in = null;
		
		try {
			in = new FileInputStream(file); 
			addFromStream(configuration, in, requiredVersion);
			
		} finally {
			try {
				in.close();
			} catch (Exception e) {
				// ignore
			}
		}
		
	}
	
	public static void addFromStream(ConfigurationModule configuration, InputStream stream, int requiredVersion) throws SAXException, IOException, ParserConfigurationException, OldConfigurationFormatException {
		Document document = XmlUtil.getInstance().parseReader(new InputStreamReader(stream));
		addFromXml(configuration, document, requiredVersion);
	}
	
	public static void addFromXml(ConfigurationModule configuration, Document xml, int requiredVersion) throws OldConfigurationFormatException {
		// check version
		if (requiredVersion > 0) {
			String contentVersion = "";
			if (xml.getDocumentElement().hasAttribute(VERSION_ATTRIBUTE)) {
				contentVersion = xml.getDocumentElement().getAttribute(VERSION_ATTRIBUTE);
			}
			if (contentVersion.equals("")) {
				throw new OldConfigurationFormatException("version not specified");
			}
			int version;
			try {
				version = Integer.parseInt(contentVersion);
			} catch (NumberFormatException e) {
				// ignore
				version = -1;
			}
			if (version < requiredVersion) {
				throw new OldConfigurationFormatException("too old version (" + version + " when required is " + requiredVersion + ")");
			}
		}

		// handle content
		addFromNode(configuration, xml.getDocumentElement());
	}
	
	private static void addFromNode(ConfigurationModule configuration, Element element) {
				
		// copy values
		NodeList entries = element.getElementsByTagName("entry");
		for (int i = 0; i < entries.getLength(); i++) {
			String name = ((Element)entries.item(i)).getAttribute("entryKey");
			NodeList values = ((Element)entries.item(i)).getElementsByTagName("value");
			for (int v = 0; v < values.getLength(); v++) {
				configuration.addValue(name, values.item(v).getTextContent());
			}
		}
		
		// recurse into modules
		NodeList modules = element.getElementsByTagName("configuration-module");
		for (int i = 0; i < modules.getLength(); i++) {
			String name = ((Element)modules.item(i)).getAttribute("moduleId");
			if (!configuration.hasSubModule(name)) {
				configuration.createSubModule(name);
			}
			addFromNode(configuration.getModule(name), (Element)modules.item(i));
		}		
	}
}
