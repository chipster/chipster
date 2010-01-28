package fi.csc.microarray.proto.repository.dummy;

import java.io.IOException;
import java.io.InputStream;
import java.util.LinkedList;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.w3c.dom.Text;
import org.xml.sax.SAXException;

import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.proto.repository.Experiment;
import fi.csc.microarray.proto.repository.Query;
import fi.csc.microarray.proto.repository.schema.ParameterClass;
import fi.csc.microarray.proto.repository.schema.ParameterInstance;
import fi.csc.microarray.proto.repository.RepositoryBase;

public class DummyRepository extends RepositoryBase {

	private ParameterClass rootClass;
	private int classHierarchyDepth;
	
	@Override
	public String getType() {
		return "Dummy repository";
	}

	@Override
	public String getIdentifier() {
		return ""; // nothing to identify, it's always the same
	}

	@Override
	public ParameterClass getRootClass() throws MicroarrayException {
		if (rootClass == null) {
			try {
				DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
				DocumentBuilder builder = factory.newDocumentBuilder();
				InputStream is = DummyRepository.class.getResourceAsStream("/MOSubClasses.xml");
				Document document = builder.parse(is);
				
				Node rootNode = document.getElementsByTagName("Classes").item(0);
				if (rootNode instanceof Element) {
					rootClass = createNamedParameterClassFromElement(
							"Root", (Element) rootNode);
					classHierarchyDepth = rootClass.setGrade(0)+1;
				} else {
					throw new IllegalArgumentException("rootNode is " + rootNode.getClass().getSimpleName() + ", value " + rootNode.toString());
				}

			} catch (ParserConfigurationException e) {
				throw new MicroarrayException(e);
			} catch (SAXException e) {
				throw new MicroarrayException(e);
			} catch (IOException e) {
				throw new MicroarrayException(e);
			}
		}
		return rootClass;
	}
	
	@Override
	public int getClassHierarchyDepth() {
		return classHierarchyDepth;
	}

	@Override
	public Iterable<Experiment> executeSubQuery(Query query) {
		DummyExperiment exp = new DummyExperiment();
		LinkedList<Experiment> exps = new LinkedList<Experiment>();
		exps.add(exp);		
		return exps;
	}
	
	@Override
	public String toString() {
		return this.getType();
	}
	
	//
	// Static machinery for parsing the XML:
	//
	
	protected static ParameterClass createParameterClassFromElement(Element classElement) {
		String name = classElement.getAttributeNS(null, "name");
		return createNamedParameterClassFromElement(name, classElement);
	}

	protected static ParameterClass createNamedParameterClassFromElement(String name, Element classElement) {	
		ParameterClass paramClass = new ParameterClass(name);
		
		NodeList childNodeList = classElement.getChildNodes();

		for (int i = 0; i < childNodeList.getLength(); i++) {
			Node childNode = childNodeList.item(i);
			
			if (childNode instanceof Element) {
				Element childElement = (Element) childNode;
				String tag = childElement.getTagName();
				
				if (tag.equals("Description")) {
					NodeList descriptionChildList = childElement.getChildNodes();
					for (int j = 0; j < descriptionChildList.getLength(); j++) {
						Node descriptionChild = descriptionChildList.item(j);
						if (descriptionChild instanceof Text) {
							Text descriptionText = (Text) descriptionChild;
							String description = parseTextElement(descriptionText);
							paramClass.setDescription(description);
						}
					}
				}
				
				else if (tag.equals("Class")) {
					ParameterClass paramSubClass = createParameterClassFromElement(childElement);
					if (paramSubClass != null) {
						paramClass.addSubclass(paramSubClass);
					}
				}
				
				else if (tag.equals("Instances")) {
					NodeList instanceList = childElement.getElementsByTagName("Instance");
					for (int j = 0; j < instanceList.getLength(); j++) {
						Element instanceElement = (Element) instanceList.item(j);
						ParameterInstance paramInstance = createParameterInstanceFromElement(instanceElement);
						paramClass.addInstance(paramInstance);
					}
				}
			}
		}
		
		return paramClass;
	}
	
	protected static ParameterInstance createParameterInstanceFromElement(Element instanceElement) {
		String name = instanceElement.getAttributeNS(null, "name");
		
		ParameterInstance paramInstance = new ParameterInstance(name);
		
		NodeList synonymList = instanceElement.getElementsByTagName("Synonym");
		for (int i = 0; i < synonymList.getLength(); i++) {
			Element synonymElement = (Element) synonymList.item(i);
			String synonymName = synonymElement.getAttributeNS(null, "name");
			paramInstance.addSynonym(synonymName);
		}
		
		return paramInstance;
	}
	
	protected static String parseTextElement(Text textElement) {
		/*
		 * Get a textual representation of the org.w3c.dom.Text element.
		 * (Strangely enough, the interface's own getWholeText()
		 * method wouldn't work, but throw AbstractMethodException.)
		 */
		String text = textElement.toString();

		/*
		 * Remove leading "[#text: " and trailing "]" substrings, which
		 * are always present in the Text object's toString value.
		 */
		text = text.substring(8, text.length()-1);
		
		/*
		 * Remove all line breaks.
		 */
		if (text.indexOf('\n') != -1) {
			text = text.replace('\n', ' ');
		}
		
		/*
		 * Finally, remove all redundant whitespace.
		 */
		while (text.indexOf("  ") != -1) {
			text = text.replace("  ", " ");
		}
		
		return text;
	}
}
