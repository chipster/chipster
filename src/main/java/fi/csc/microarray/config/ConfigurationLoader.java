package fi.csc.microarray.config;

import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import fi.csc.microarray.util.XmlUtil;

public class ConfigurationLoader {
	
	private static final String VERSION_ATTRIBUTE = "content-version";
	
	public static class IllegalConfigurationException extends Exception {

		public IllegalConfigurationException(String string) {
			super(string);
		}
		
	}
	
	private int requiredVersion;
	private Configuration configuration;
	
	public ConfigurationLoader(Configuration configuration, int requiredVersion) {
		
		this.configuration = configuration;
		this.requiredVersion = requiredVersion;
	}

	public void addFromStream(InputStream stream, boolean isSpecification) throws SAXException, IOException, ParserConfigurationException, IllegalConfigurationException {
		Document document = XmlUtil.parseReader(new InputStreamReader(stream));
		addFromXml(document, isSpecification);
	}
	
	public void addFromXml(Document xml, boolean isSpecification) throws IllegalConfigurationException {
		// check version
		if (requiredVersion > 0) {
			String contentVersion = "";
			if (xml.getDocumentElement().hasAttribute(VERSION_ATTRIBUTE)) {
				contentVersion = xml.getDocumentElement().getAttribute(VERSION_ATTRIBUTE);
			}
			if (contentVersion.equals("")) {
				throw new IllegalConfigurationException("version not specified");
			}
			int version;
			try {
				version = Integer.parseInt(contentVersion);
			} catch (NumberFormatException e) {
				// ignore
				version = -1;
			}
			if (version < requiredVersion) {
				throw new IllegalConfigurationException("too old version (" + version + " when required is " + requiredVersion + ")");
			}
		}

		// handle content
		addFromNode(configuration.getRootModule(), xml.getDocumentElement(), isSpecification);
	}
	
	private void addFromNode(ConfigurationModule module, Element element, boolean isSpecification) throws IllegalConfigurationException {
				
		// copy values
		NodeList entries = element.getChildNodes();
		for (int i = 0; i < entries.getLength(); i++) {
			Node item = entries.item(i);
			if (item instanceof Element && "entry".equals(item.getNodeName())) {
				Element entry = (Element)item;
				String name = entry.getAttribute("entryKey");
				String type = entry.getAttribute("type");
				NodeList values = entry.getElementsByTagName("value");

				if (name == null || type == null) {
					throw new IllegalConfigurationException("missing entryKey or type for " + name);
				}

				// locate our entry, possibly creating a new one
				ConfigurationEntry configurationEntry;
				if (isSpecification) {
					configurationEntry = new ConfigurationEntry(name, type);
					module.addEntry(configurationEntry);
				} else {
					configurationEntry = module.getEntry(name);
					if (configurationEntry == null) {
						throw new IllegalConfigurationException("unsupported entry: " + name);
					}
				}
				
				// mark to be set or fill with values
				if (entry.hasAttribute("mustBeSet") && entry.getAttribute("mustBeSet").equals("true")) {
					if (values.getLength() > 0) {
						throw new IllegalConfigurationException("illegal config specification: both values and mustBeSet given for " + name);
					}
					
					configurationEntry.setMustBeSet(true);

				} else {

					String[] textValues = new String[values.getLength()];
					for (int v = 0; v < values.getLength(); v++) {
						textValues[v] = values.item(v).getTextContent();
					}
					configurationEntry.setValue(textValues);
					
				}
			}
		}
		
		// recurse into modules
		NodeList modules = element.getElementsByTagName("configuration-module");
		for (int i = 0; i < modules.getLength(); i++) {
			Element submoduleElement = (Element)modules.item(i);
			String name = submoduleElement.getAttribute("moduleId");
			if (!configuration.isModuleEnabled(name)) {
				if (!isSpecification) {
					// tried to introduce extra module => error
					throw new IllegalConfigurationException("unsupported module: " + name);
				}

				continue; // skip disabled modules
			}
			if (!module.hasSubModule(name)) {
				module.createSubModule(name);
			}
			addFromNode(module.getModule(name), submoduleElement, isSpecification);
		}		
	}
}
