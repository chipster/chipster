package fi.csc.microarray.config;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import java.net.Inet4Address;
import java.net.InetAddress;
import java.net.NetworkInterface;
import java.net.SocketException;
import java.net.UnknownHostException;
import java.util.Collections;
import java.util.Enumeration;
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

	static final String CURRENT_R_VERSION = "R-2.12";
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
			{"admin e-mail address", "chipster-admin@mydomain"},
			{"path to R binary", "/opt/chipster/tools/R/bin/R"},
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
	private final int MANAGER_EMAIL = 8;
	private final int R_COMMAND_INDEX = 9;
	private final int MAX_JOBS_INDEX = 10;

	private String[][] passwords = new String[][] {
			{"comp", ""},
			{"auth", ""},
			{"filebroker", ""},
			{"manager", ""}
	};
	
	private HashMap<String, Document> documentsToWrite = new HashMap<String, Document>();

	public ConfigTool() throws ParserConfigurationException {
		System.out.println("Configuring Chipster");
	}
	
	public static void main(String[] args) throws Exception {
		ConfigTool configTool = new ConfigTool();
		UpgradeTool upgradeTool = new UpgradeTool();
		SetupTool setupTool = new SetupTool(); 
		
		if (args.length == 0) {
			fail();

		} else if ("configure".equals(args[0])) {
			configTool.configure();
		} else if ("auto-configure".equals(args[0])) {
			configTool.simpleConfigure(null);
		} else if ("simple-configure".equals(args[0]) && args.length == 2) {
			configTool.simpleConfigure(args[1]);
		} else if ("genpasswd".equals(args[0])) {
			configTool.genpasswd();

		} else if ("setup".equals(args[0])) {
			setupTool.setup();
				
		} else if (args[0].startsWith("upgrade")) {
			String[] parts = args[0].split("_");
			int fromMajor = Integer.parseInt(parts[1]);
			int toMajor = Integer.parseInt(parts[2]);
			if (args.length > 1) {
				upgradeTool.upgrade(new File(args[1]), fromMajor, toMajor);
			} else {
				System.out.println("Please specify location of the old installation directory as an argument (e.g., \"./upgrade.sh /opt/chipster-1.2.3\")");
			}

		} else {
			fail();
		}
	}
	
	private static void fail() {
		System.out.println("Illegal arguments! Please specify one of: configure, genpasswd, setup, upgrade_<major version number of source>_<major version number of target>");
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
			XmlUtil.printXml(documentsToWrite.get(file), new OutputStreamWriter(new FileOutputStream(file)));
		}
		System.out.println("\nAll changes successfully written!");
	}

	public static void verifyChanges(BufferedReader in) throws Exception {
		verifyChanges(in, "Please verify changes. Should changes be written to disk");	
	}

	public static void verifyChanges(BufferedReader in, String question) throws Exception {
		System.out.println(question + " [yes/no]?");
		for (String answer = in.readLine();!"yes".equals(answer); answer = in.readLine()) {
			if ("no".equals(answer)) {
				throw new Exception("User decided to abort");	
			}
			System.out.println(question + " [yes/no]?");
		}
	}

	public void configure() throws Exception {

		try {
			System.out.println("No changes are written before you verify them.");
			System.out.println("");

			BufferedReader in = new BufferedReader(new InputStreamReader(System.in));

			//
			// STEP 1. GATHER DATA
			//
			
			// sniff current host
			try {
				String host = getInetAddress().getHostName();
				configs[BROKER_HOST_INDEX][VAL_INDEX] = host;
				configs[FILEBROKER_HOST_INDEX][VAL_INDEX] = host;
				configs[WS_CODEBASE_INDEX][VAL_INDEX] = "http://" + configs[BROKER_HOST_INDEX][VAL_INDEX] + ":" + configs[WS_PORT][VAL_INDEX];
			} catch (UnknownHostException e) {
				// ignore, sniffing failed
			}
			
			// gather required data
			for (int i = 0; i < configs.length; i++) {
				if (i == WS_CODEBASE_INDEX) {
					continue;
				}
				System.out.println("Please specify " + configs[i][KEY_INDEX] + " [" + configs[i][VAL_INDEX] + "]: ");
				String line = in.readLine();
				if (!line.trim().equals("")) {
					configs[i][VAL_INDEX] = line;
				}
			}
			
			// add web start location
			configs[WS_CODEBASE_INDEX][VAL_INDEX] = "http://" + configs[BROKER_HOST_INDEX][VAL_INDEX] + ":" + configs[WS_PORT][VAL_INDEX];

			
			
			
			//
			// STEP 2. UPDATE CONFIGS
			//
			updateConfigs();
			
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

	/**
	 * Configure server host, use defaults for other things.
	 * 
	 * @param host
	 * @throws Exception
	 */
	public void simpleConfigure(String host) throws Exception {
	
		// auto detect hostname
		if (host == null) {
			host = getInetAddress().getHostName();
		}
		
		try {

			configs[BROKER_HOST_INDEX][VAL_INDEX] = host;
			configs[FILEBROKER_HOST_INDEX][VAL_INDEX] = host;
			configs[WS_CODEBASE_INDEX][VAL_INDEX] = "http://" + host + ":" + configs[WS_PORT][VAL_INDEX];
			
			updateConfigs();

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

	
	
	private void updateConfigs() throws Exception {
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
		File runtimesConfigFile = new File("comp" + File.separator + DirectoryLayout.CONF_DIR + File.separator + "runtimes.xml");
		if (runtimesConfigFile.exists()) {
			updateRuntimesConfigFile(runtimesConfigFile);
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
	}

	private void updateActivemqConfigFile(File configFile) throws Exception {
		Document doc = openForUpdating("ActiveMQ", configFile);
		Element broker = (Element)doc.getDocumentElement().getElementsByTagName("broker").item(0);
		
		Element transportConnectors = (Element)broker.getElementsByTagName("transportConnectors").item(0);		
		Element transportConnector = (Element)transportConnectors.getElementsByTagName("transportConnector").item(0); // edit first in the list (could use attribute name to decide right one)..
		String uri = configs[BROKER_PROTOCOL_INDEX][VAL_INDEX] + "://" + configs[BROKER_HOST_INDEX][VAL_INDEX] + ":" + configs[BROKER_PORT_INDEX][VAL_INDEX];
		updateElementAttribute(transportConnector, "uri", uri);
		
		writeLater(configFile, doc);
	}
	
	private void updateWsConfigFile(File configFile) throws Exception {
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

	private void updateActivemqConfigFilePasswords(File configFile) throws Exception {
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

		Element securityModule = XmlUtil.getChildWithAttributeValue(doc.getDocumentElement(), "moduleId", "security");
		Element usernameElement = XmlUtil.getChildWithAttributeValue(securityModule, "entryKey", "username");
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

		Element messagingModule = XmlUtil.getChildWithAttributeValue(doc.getDocumentElement(), "moduleId", "messaging");
		updateConfigEntryValue(messagingModule, "broker-host", configs[BROKER_HOST_INDEX][VAL_INDEX]);
		updateConfigEntryValue(messagingModule, "broker-protocol", configs[BROKER_PROTOCOL_INDEX][VAL_INDEX]);
		updateConfigEntryValue(messagingModule, "broker-port", configs[BROKER_PORT_INDEX][VAL_INDEX]);

		Element filebrokerModule = XmlUtil.getChildWithAttributeValue(doc.getDocumentElement(), "moduleId", "filebroker");
		if (filebrokerModule != null) {
			updateConfigEntryValue(filebrokerModule, "port", configs[FILEBROKER_PORT_INDEX][VAL_INDEX]);
			updateConfigEntryValue(filebrokerModule, "url", createFilebrokerUrl());
		}

		Element analyserModule = XmlUtil.getChildWithAttributeValue(doc.getDocumentElement(), "moduleId", "comp");
		if (analyserModule != null) {
			updateConfigEntryValue(analyserModule, "max-jobs", configs[MAX_JOBS_INDEX][VAL_INDEX]);
		}
		
		Element webstartModule = XmlUtil.getChildWithAttributeValue(doc.getDocumentElement(), "moduleId", "webstart");
		if (webstartModule != null) {
			updateConfigEntryValue(webstartModule, "port", configs[WS_PORT][VAL_INDEX]);
		}

		Element clientModule = XmlUtil.getChildWithAttributeValue(doc.getDocumentElement(), "moduleId", "client");
		if (clientModule != null) {
			updateConfigEntryValue(clientModule, "manual-root", configs[WS_CODEBASE_INDEX][VAL_INDEX] + "/manual/");
		}

		
		Element managerModule = XmlUtil.getChildWithAttributeValue(doc.getDocumentElement(), "moduleId", "manager");
		if (managerModule != null) {
			updateConfigEntryValue(managerModule, "web-console-port", configs[MANAGER_PORT][VAL_INDEX]);
			updateConfigEntryValue(managerModule, "admin-email", configs[MANAGER_EMAIL][VAL_INDEX]);
		}

		writeLater(configFile, doc);
	}

	private void updateRuntimesConfigFile(File configFile) throws Exception {
		
		boolean ok = false;
		Document doc = openForUpdating("Runtimes", configFile);
		Element runtimesElement = (Element)doc.getElementsByTagName("runtimes").item(0);
		for (Element runtimeElement: XmlUtil.getChildElements(runtimesElement, "runtime")) {
			String runtimeName = XmlUtil.getChildElement(runtimeElement, "name").getTextContent();
			if (runtimeName.equals(CURRENT_R_VERSION)) {
				Element handlerElement = XmlUtil.getChildElement(runtimeElement, "handler");
				for (Element parameterElement: XmlUtil.getChildElements(handlerElement, "parameter")) {
					String paramName = XmlUtil.getChildElement(parameterElement, "name").getTextContent();
					if (paramName.equals("command")) {
						Element commandValueElement = XmlUtil.getChildElement(parameterElement, "value");
						updateElementValue(commandValueElement, "R-2.6.1 command", configs[R_COMMAND_INDEX][VAL_INDEX]);
						ok = true;
					}
				}
			} 
		}

		if (ok) {
			writeLater(configFile, doc);
		} else {
			throw new RuntimeException("Could not update R-2.6.1 command to runtimes.xml");
		}
	}

	
	
	private String createFilebrokerUrl() {
		return "http://" + configs[FILEBROKER_HOST_INDEX][VAL_INDEX];
	}

	private void updateConfigEntryValue(Element module, String name, String newValue) {
		Element entry = XmlUtil.getChildWithAttributeValue(module, "entryKey", name);
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

	private Document openForUpdating(String name, File configFile) throws SAXException, IOException, ParserConfigurationException {
		System.out.println("Updating " + name + " config in " + configFile.getAbsolutePath());
		Document doc = XmlUtil.parseFile(configFile);
		return doc;
	}

	public static String[] getComponentDirsWithConfig() {
		return componentDirsWithConfig;
	}

	
	/**
	 * Prefer something other than loopback.
	 * 
	 * @return
	 * @throws SocketException
	 * @throws UnknownHostException
	 */
	public static InetAddress getInetAddress() throws SocketException, UnknownHostException {
		
		// try to use eth0 if found
		NetworkInterface eth0 = NetworkInterface.getByName("eth0");
		if (eth0 != null) {
			Enumeration<InetAddress> addresses = eth0.getInetAddresses();
			for (InetAddress address : Collections.list(addresses)) {
				if (address instanceof Inet4Address) {
					return address;
				}
			}
		} 
		
		// something other than eth0 or loopback
		Enumeration<NetworkInterface> interfaces = NetworkInterface.getNetworkInterfaces();
		for (NetworkInterface iface : Collections.list(interfaces)) {
			Enumeration<InetAddress> addresses = iface.getInetAddresses();
			for (InetAddress address : Collections.list(addresses)) {
				if (address instanceof Inet4Address && !address.isLoopbackAddress()) {
					return address;
				}
			}
		}
	
		// default, probably loopback
		return InetAddress.getLocalHost();
	}

}
