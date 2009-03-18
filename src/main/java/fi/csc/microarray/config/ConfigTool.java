package fi.csc.microarray.config;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import java.net.InetAddress;
import java.net.UnknownHostException;
import java.util.HashMap;
import java.util.UUID;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.TransformerException;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import fi.csc.microarray.util.XmlUtil;

/**
 * Simple tool for centrally changing configuration of the Chipster server environment.
 * 
 * @author Aleksi Kallio
 *
 */
public class ConfigTool {

	private final String brokerDir = "activemq";
	private final String webstartDir = "webstart";

	private final static String[] componentDirsWithConfig = new String[] {
			"comp",
			"auth",
			"fileserver",
			"manager",
			"client",
			"webstart"
	};

	private String[][] configs = new String[][] {
			{"message broker (ActiveMQ) host", "myhost.mydomain"},
			{"message broker protocol", "tcp"},
			{"message broker port", "61616"},
			{"file broker host", "myhost.mydomain"},
			{"file broker port", "8080"},
			{"URL of Web Start files", "http://myhost.mydomain"},
			{"Web Start www-server port", "8081"},
			{"manager www-console port", "8082"},
			{"R command", "R"},
			{"max. simultanous jobs (more recommended when compute service on separate node)", "3"}
	};

	private final int KEY_INDEX = 0;
	private final int VAL_INDEX = 1;

	private final int BROKER_HOST_INDEX = 0;
	private final int BROKER_PROTOCOL_INDEX = 1;
	private final int BROKER_PORT_INDEX = 2;
	private final int FILEBROKER_HOST_INDEX = 3;
	private final int FILEBROKER_PORT_INDEX = 4;
	private final int WS_CODEBASE_INDEX = 5;
	private final int WS_PORT = 6;
	private final int MANAGER_PORT = 7;
	private final int R_COMMAND_INDEX = 8;
	private final int MAX_JOBS_INDEX = 9;

	private String[][] passwords = new String[][] {
			{"comp", ""},
			{"auth", ""},
			{"filebroker", ""},
			{"manager", ""}
	};
	
	private HashMap<String, Document> documentsToWrite = new HashMap<String, Document>();

	private XmlUtil xml;

	public ConfigTool() throws ParserConfigurationException {
		this.xml = XmlUtil.getInstance();
		System.out.println("Chipster ConfigTool");
		System.out.println("");
		System.out.println("No changes are written before you verify them");
		System.out.println("");
	}
	
	public static void main(String[] args) throws Exception {
		ConfigTool configTool = new ConfigTool();
		UpgradeTool upgradeTool = new UpgradeTool();
		
		if (args.length == 0) {
			fail();

		} else if ("configure".equals(args[0])) {
			configTool.configure();

		} else if ("genpasswd".equals(args[0])) {
			configTool.genpasswd();

		} else if ("upgrade".equals(args[0])) {
			if (args.length > 1) {
				upgradeTool.upgrade(new File(args[1]));
			} else {
				System.out.println("Please specify location of the old installation directory as an argument (e.g., \"./upgrade.sh /opt/chipster-1.2.3\")");
			}

		} else {
			fail();
		}
	}
	
	private static void fail() {
		System.out.println("Illegal arguments! Please specify one of: configure, genpasswd, upgrade");
	}
	
	private void genpasswd() throws Exception {

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(System.in));

			//
			// STEP 1. GATHER DATA
			//

			// generate passwords
			for (int i = 0; i < passwords.length; i++) {
				passwords[i][VAL_INDEX] = UUID.randomUUID().toString(); 
			}
			
			//
			// STEP 2. UPDATE CONFIGS
			//
			
			// update all Chipster configs
			for (String componentDir : getComponentDirsWithConfig()) {
				if (new File(componentDir).exists()) {
					File configFile = new File(componentDir + File.separator + DirectoryLayout.CONF_DIR + File.separator + Configuration.CONFIG_FILENAME);
					updateChipsterConfigFilePasswords(configFile);
				}
			}

			// update ActiveMQ config
			File activemqConfigFile = new File(brokerDir + File.separator + DirectoryLayout.CONF_DIR + File.separator + "activemq.xml");
			if (activemqConfigFile.exists()) {
				updateActivemqConfigFilePasswords(activemqConfigFile);
			}
			
			verifyChanges(in);

		} catch (Throwable t) {
			t.printStackTrace();
			System.err.println("\nQuitting, no changes written to disk!");
			return;

		}
		
		
		//
		// STEP 3. WRITE CHANGES
		//
		
		writeChangesToDisk();
		
	}

	private void writeChangesToDisk() throws TransformerException, UnsupportedEncodingException, FileNotFoundException {
		// write out files
		for (String file : documentsToWrite.keySet()) {
			System.out.println("Writing changes to " + file + "...");
			xml.printXml(documentsToWrite.get(file), new OutputStreamWriter(new FileOutputStream(file)));
		}
		System.out.println("\nAll changes successfully written!");
	}

	public static void verifyChanges(BufferedReader in) throws Exception {
		System.out.println("Please verify changes. Should changes be written to disk [yes/no]?");
		String answer = in.readLine();
		if (!"yes".equals(answer)) {
			throw new Exception("User decided to abort");
		}
	}

	public void configure() throws Exception {

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(System.in));

			//
			// STEP 1. GATHER DATA
			//
			
			// sniff current host
			try {
				String host = InetAddress.getLocalHost().getHostName();
				configs[BROKER_HOST_INDEX][VAL_INDEX] = host;
				configs[FILEBROKER_HOST_INDEX][VAL_INDEX] = host;
				configs[WS_CODEBASE_INDEX][VAL_INDEX] = "http://" + host + ":8081";
			} catch (UnknownHostException e) {
				// ignore, sniffing failed
			}
			
			// gather required data
			for (int i = 0; i < configs.length; i++) {
				System.out.println("Please specify " + configs[i][KEY_INDEX] + " [" + configs[i][VAL_INDEX] + "]: ");
				String line = in.readLine();
				if (!line.trim().equals("")) {
					configs[i][VAL_INDEX] = line;
				}
			}

			//
			// STEP 2. UPDATE CONFIGS
			//
			
			// update all Chipster configs
			for (String componentDir : getComponentDirsWithConfig()) {
				if (new File(componentDir).exists()) {
					File configFile = new File(componentDir + File.separator + DirectoryLayout.CONF_DIR + File.separator + Configuration.CONFIG_FILENAME);
					updateChipsterConfigFile(configFile);
				}
			}
			File wsClientConfigFile = new File("webstart" + File.separator + DirectoryLayout.WEB_ROOT + File.separator + Configuration.CONFIG_FILENAME);
			if (wsClientConfigFile.exists()) {
				updateChipsterConfigFile(wsClientConfigFile);
			}

			// update ActiveMQ config
			File activemqConfigFile = new File(brokerDir + File.separator + DirectoryLayout.CONF_DIR + File.separator + "activemq.xml");
			if (activemqConfigFile.exists()) {
				updateActivemqConfigFile(activemqConfigFile);
			}
			
			// update Web Start config
			File wsConfigFile = new File(webstartDir + File.separator + DirectoryLayout.WEB_ROOT + File.separator + "chipster.jnlp");
			if (wsConfigFile.exists()) {
				updateWsConfigFile(wsConfigFile);
			}
			
			verifyChanges(in);

		} catch (Throwable t) {
			t.printStackTrace();
			System.err.println("\nQuitting, no changes written to disk!");
			return;

		}
		
		
		//
		// STEP 3. WRITE CHANGES
		//
		
		writeChangesToDisk();

	}

	private void updateWsConfigFile(File configFile) throws SAXException, IOException, TransformerException, UnsupportedEncodingException, FileNotFoundException {
		Document doc = openForUpdating("Web Start", configFile);
		Element jnlp = (Element)doc.getDocumentElement();
		updateElementAttribute(jnlp, "codebase", configs[WS_CODEBASE_INDEX][VAL_INDEX]);
		Element applicationDesc = (Element)jnlp.getElementsByTagName("application-desc").item(0);
		NodeList arguments = applicationDesc.getElementsByTagName("argument");
		Element lastArgument = (Element)arguments.item(arguments.getLength() - 1);
		String url = "http://" + configs[BROKER_HOST_INDEX][VAL_INDEX] + ":" + configs[WS_PORT][VAL_INDEX] + "/" + Configuration.CONFIG_FILENAME;
		updateElementValue(lastArgument, "configuration URL (for Web Start)", url);
		writeLater(configFile, doc);
	}

	private void updateActivemqConfigFile(File configFile) throws SAXException, IOException, TransformerException, UnsupportedEncodingException, FileNotFoundException {
		Document doc = openForUpdating("ActiveMQ", configFile);
		Element broker = (Element)doc.getDocumentElement().getElementsByTagName("broker").item(0);
		
		Element transportConnectors = (Element)broker.getElementsByTagName("transportConnectors").item(0);		
		Element transportConnector = (Element)transportConnectors.getElementsByTagName("transportConnector").item(0); // edit first in the list (could use attribute name to decide right one)..
		String uri = configs[BROKER_PROTOCOL_INDEX][VAL_INDEX] + "://" + configs[BROKER_HOST_INDEX][VAL_INDEX] + ":" + configs[BROKER_PORT_INDEX][VAL_INDEX];
		updateElementAttribute(transportConnector, "uri", uri);
		
		writeLater(configFile, doc);
	}
	
	private void updateActivemqConfigFilePasswords(File configFile) throws SAXException, IOException, TransformerException, UnsupportedEncodingException, FileNotFoundException {
		Document doc = openForUpdating("ActiveMQ", configFile);
		Element broker = (Element)doc.getDocumentElement().getElementsByTagName("broker").item(0);
			
		NodeList users = ((Element)((Element)((Element)broker.getElementsByTagName("plugins").item(0)).getElementsByTagName("simpleAuthenticationPlugin").item(0)).getElementsByTagName("users").item(0)).getElementsByTagName("authenticationUser");
		for (int i = 0; i < users.getLength(); i++) {
			for (int p = 0; p < passwords.length; p++) {
				Element user = (Element)users.item(i);
				if (user.getAttribute("username").equals(passwords[p][KEY_INDEX])) {
					updateElementAttribute(user, "password for " + passwords[p][KEY_INDEX], "password", passwords[p][VAL_INDEX]);
					break;
				}
			}
		}
		writeLater(configFile, doc);
	
	}

	private void updateChipsterConfigFilePasswords(File configFile) throws Exception {
		Document doc = openForUpdating("Chipster", configFile);

		Element securityModule = xml.getChildWithAttribute(doc.getDocumentElement(), "moduleId", "security");
		Element usernameElement = xml.getChildWithAttribute(securityModule, "entryKey", "username");
		String username = ((Element)usernameElement.getElementsByTagName("value").item(0)).getTextContent();
		for (int i = 0; i < passwords.length; i++) {
			if (username.equals(passwords[i][KEY_INDEX])) {
				updateConfigEntryValue(securityModule, "password", passwords[i][VAL_INDEX]);
				break;
			}
		}
		writeLater(configFile, doc);
	}
	
	private void updateChipsterConfigFile(File configFile) throws Exception {
		Document doc = openForUpdating("Chipster", configFile);

		Element messagingModule = xml.getChildWithAttribute(doc.getDocumentElement(), "moduleId", "messaging");
		updateConfigEntryValue(messagingModule, "broker-host", configs[BROKER_HOST_INDEX][VAL_INDEX]);
		updateConfigEntryValue(messagingModule, "broker-protocol", configs[BROKER_PROTOCOL_INDEX][VAL_INDEX]);
		updateConfigEntryValue(messagingModule, "broker-port", configs[BROKER_PORT_INDEX][VAL_INDEX]);

		Element filebrokerModule = xml.getChildWithAttribute(doc.getDocumentElement(), "moduleId", "filebroker");
		if (filebrokerModule != null) {
			updateConfigEntryValue(filebrokerModule, "port", configs[FILEBROKER_PORT_INDEX][VAL_INDEX]);
			updateConfigEntryValue(filebrokerModule, "url", createFilebrokerUrl());
		}

		Element analyserModule = xml.getChildWithAttribute(doc.getDocumentElement(), "moduleId", "comp");
		if (analyserModule != null) {
			updateConfigEntryValue(analyserModule, "max-jobs", configs[MAX_JOBS_INDEX][VAL_INDEX]);
			updateConfigEntryValue(analyserModule, "r-command", configs[R_COMMAND_INDEX][VAL_INDEX]);
		}
		
		Element webstartModule = xml.getChildWithAttribute(doc.getDocumentElement(), "moduleId", "webstart");
		if (webstartModule != null) {
			updateConfigEntryValue(webstartModule, "port", configs[WS_PORT][VAL_INDEX]);
		}

		Element managerModule = xml.getChildWithAttribute(doc.getDocumentElement(), "moduleId", "manager");
		if (managerModule != null) {
			updateConfigEntryValue(managerModule, "web-console-port", configs[MANAGER_PORT][VAL_INDEX]);
		}

		writeLater(configFile, doc);
	}

	private String createFilebrokerUrl() {
		return "http://" + configs[FILEBROKER_HOST_INDEX][VAL_INDEX];
	}

	private void updateConfigEntryValue(Element module, String name, String newValue) {
		Element entry = xml.getChildWithAttribute(module, "entryKey", name);
		Element value = (Element)entry.getElementsByTagName("value").item(0);
		updateElementValue(value, name, newValue);
	}

	private void updateElementValue(Element element, String logicalName, String newValue) {
		System.out.println("  changing " + logicalName + ": " + element.getTextContent() + " -> " + newValue);
		element.setTextContent(newValue);
	}

	private void updateElementAttribute(Element element, String attrName, String attrValue) {
		updateElementAttribute(element, attrName, attrName, attrValue);
	}
	
	private void updateElementAttribute(Element element, String logicalName, String attrName, String attrValue) {
		System.out.println("  changing " + logicalName + ": " + element.getAttribute(attrName) + " -> " + attrValue);
		element.setAttribute(attrName, attrValue);
	}

	private void writeLater(File configFile, Document doc) throws TransformerException, UnsupportedEncodingException, FileNotFoundException {
		documentsToWrite.put(configFile.getAbsolutePath(), doc);
		System.out.println("");
	}

	private Document openForUpdating(String name, File configFile) throws SAXException, IOException {
		System.out.println("Updating " + name + " config in " + configFile.getAbsolutePath());
		Document doc = xml.parseFile(configFile);
		return doc;
	}

	public static String[] getComponentDirsWithConfig() {
		return componentDirsWithConfig;
	}
	
}
