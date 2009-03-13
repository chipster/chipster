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
		Document document = XmlUtil.getInstance().parseReader(new InputStreamReader(stream));
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
				NodeList values = entry.getElementsByTagName("value");

				if (entry.getElementsByTagName("mustBeSet").getLength() > 0) {
					if (values.getLength() > 0) {
						throw new RuntimeException("illegal config specification: both value and mustBeSet given for " + name);
					}
					module.putValues(name, ConfigurationModule.VALUE_MUST_BE_SET);

				} else {

					
					if (module.getValues(name) != null && module.getValues(name) != ConfigurationModule.VALUE_MUST_BE_SET && !isSpecification) {
						System.out.println("Warning: overriding " + name);
					}
					
					if (module.getValues(name) != null) {
						// replace old valueset, possibly the special value ConfigurationModule.VALUE_MUST_BE_SET
						module.removeValues(name);
					} else if (!isSpecification) {
						// is not a replace and new values are not allowed => error
						throw new IllegalConfigurationException("unsupported entry: " + name);
					}

					for (int v = 0; v < values.getLength(); v++) {
						module.addValue(name, values.item(v).getTextContent());
					}
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
