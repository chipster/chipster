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
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.List;
import java.util.UUID;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.TransformerException;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import com.sun.org.apache.xerces.internal.util.URI;

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
			"webstart",
			"jobmanager",
			"toolbox"
	};

	private String[][] configs = new String[][] {
			{"broker public host/ip", "myhost.mydomain"},
			{"broker private host/ip", "myhost.mydomain"},
			{"message broker protocol", "ssl"},
			{"message broker port", "61616"},
			{"file broker protocol", "https"},
			{"file broker port", "8080"},
			{"Web Start www-server port", "8081"},
			{"manager www-console port", "8082"},
			{"admin e-mail address", "chipster-admin@mydomain"},
			{"max. simultanous jobs (more recommended when compute service on separate node)", "3"}
	};

	private final int KEY_INDEX = 0;
	private final int VAL_INDEX = 1;

	private final int BROKER_PUBLIC_HOST_INDEX = 0;
	private final int BROKER_PRIVATE_HOST_INDEX = 1;
	private final int BROKER_PROTOCOL_INDEX = 2;
	private final int BROKER_PORT_INDEX = 3;
	private final int FILEBROKER_PORT_PROTOCOL = 4;
	private final int FILEBROKER_PORT_INDEX = 5;
	private final int WS_PORT = 6;
	private final int MANAGER_PORT = 7;
	private final int MANAGER_EMAIL = 8;
	private final int MAX_JOBS_INDEX = 9;

	private String[][] passwords = new String[][] {
			{"comp", ""},
			{"auth", ""},
			{"filebroker", ""},
			{"manager", ""},
			{"jobmanager", ""},
			{"toolbox", ""}
	};
	
	private HashMap<String, Document> documentsToWrite = new HashMap<>();
	private HashMap<String, List<String>> filesToWrite = new HashMap<>();

	public static void main(String[] args) throws Exception {
		ConfigTool configTool = new ConfigTool();

		if (args.length == 0 || "configure".equals(args[0])) {
			System.out.println("Configuring Chipster...");
			configTool.configure();
			
		} else if ("auto-configure".equals(args[0]) || "auto".equals(args[0])) {
			System.out.println("Configuring Chipster...");
			configTool.simpleConfigure(null, null);
			
		} else if ("simple-configure".equals(args[0]) && args.length == 3) {
			System.out.println("Configuring Chipster...");
			configTool.simpleConfigure(args[1], args[2]);			
			
		} else if ("simple-configure".equals(args[0]) && args.length == 2) {
			System.out.println("Configuring Chipster...");
			configTool.simpleConfigure(args[1], null);			

		} else if ("genpasswd".equals(args[0])) {
			System.out.println("Generating Chipster server password...");
			configTool.genpasswd();
			
		} else if ("edit".equals(args[0]) && args.length >= 4) {
			
			// Check components to process
			String[] components;
			if ("servers".equals(args[1])) {
				components = componentDirsWithConfig;
			} else {
				components = new String[] { args[1] };
			}
						
			// Process the edit command
			if (!"print".equals(args[2])) {
				System.out.println("Configuring Chipster...");
			}
			configTool.editConfig(components, args[2], args[3], args.length > 4 ? args[4] : null);

		} else if ("help".equals(args[0]) || "--help".equals(args[0]) || "-h".equals(args[0])) {
			printHelp();
			
		} else {
			commandNotUnderstood();
		}						
	}

	private static void commandNotUnderstood() {
		System.out.println("Incorrect parameter syntax. Use \"help\" for more information"); 
	}
	
	private static void printHelp() {
		System.out.println("Use:\n" + 
				"\n" + 
				"  configure.sh [configure]\n" + 
				"    Start interactive configuration utility\n" + 
				"\n" + 
				"  configure.sh [auto-configure | auto]\n" + 
				"    Detect IP automatically and configure everything else with default values\n" + 
				"\n" + 
				"  configure.sh simple-configure public-ip [private-ip]\n" + 
				"    Use given IP address (and separate private IP, when given) and configure \n" + 
				"    everything else with default values\n" + 
				"\n" + 
				"  configure.sh genpasswd\n" + 
				"    Generate strong random passwords for connections between server components \n" + 
				"    and message broker (AMQ).\n" + 
				"\n" + 
				"  configure.sh edit [<component name> | all] set entry-name entry-value\n" + 
				"    The \"edit\" command makes small changes to chipster-config.xml files and it\n" + 
				"    is mostly provided as an interface to use from automatic configuration\n" + 
				"    management and installation tools.\n" + 
				"    Set the value of entry with given name to given value. If the entry does not \n" + 
				"    exist, it is created. Entry name follows pattern CATEGORY_NAME/ENTRY_NAME,\n" + 
				"    such as comp/max-jobs. Component name is one of the available components, \n" + 
				"    \"servers\" to change all server components or \"client\"\n" + 
				"\n" + 
				"  configure.sh edit [<component name> | all] remove entry-name\n" + 
				"    Remove entry with given name completely. If the entry has default value, the\n" + 
				"    system will proceed to use it. If the entry does not exist nothing is done.\n" + 
				"\n" + 
				"  configure.sh edit [<component name> | all] print entry-name\n" + 
				"    Print the value of entry with given name to STDOUT. If the entry does not\n" + 
				"    exist nothing is printed.\n" + 
				"\n" + 
				"  configure.sh [help | --help | -h]\n" + 
				"    Print this help text.\n" + 
				"");
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

	private void writeChangesToDisk() throws TransformerException, IOException {
		// write out files
		for (String file : documentsToWrite.keySet()) {
			System.out.println("Writing changes to: " + file);
			XmlUtil.printXml(documentsToWrite.get(file), new OutputStreamWriter(new FileOutputStream(file)));
		}
		
		for (String file : filesToWrite.keySet()) {
			System.out.println("Writing changes to: " + file);			
			Files.write(Paths.get(new File(file).toURI()), filesToWrite.get(file), Charset.defaultCharset());
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
				configs[BROKER_PUBLIC_HOST_INDEX][VAL_INDEX] = host;
				configs[BROKER_PRIVATE_HOST_INDEX][VAL_INDEX] = host;
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
	public void simpleConfigure(String publicHost, String privateHost) throws Exception {
	
		// auto detect hostname(s)
		if (publicHost == null) {
			publicHost = getInetAddress().getHostName();
		}
		if (privateHost == null) {
			privateHost = publicHost;
		}
		
		try {
			configs[BROKER_PUBLIC_HOST_INDEX][VAL_INDEX] = publicHost;
			configs[BROKER_PRIVATE_HOST_INDEX][VAL_INDEX] = privateHost;
			
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


	private void editConfig(String[] components, String command, String entryName, String entryValue) throws Exception {

		try {

			// Collect changes for each config file
			for (String component : components) {
				File configFile = null;
				if ("client".equals(component)) {
					configFile = getClientConfigFile();
				} else {
					configFile = new File(component + File.separator + DirectoryLayout.CONF_DIR + File.separator + Configuration.CONFIG_FILENAME);
				}
				Document doc = openXmlForUpdating("Chipster", configFile, true);

				String[] nameParts = entryName.split("/");
				if (nameParts.length != 2) {
					throw new IllegalArgumentException("illegal entry name: " + entryName);
				}
				String category = nameParts[0];
				String name = nameParts[1];
				Element module = XmlUtil.getChildWithAttributeValue(doc.getDocumentElement(), "moduleId", category);
				if (module == null) {
					throw new IllegalArgumentException("illegal category name: " + category);
				}
				if ("set".equals(command)) {
					setConfigEntryValue(doc, module, name, entryValue);
					writeLater(configFile, doc);
					
				} else if ("remove".equals(command)) {
					removeConfigEntry(module, name);
					writeLater(configFile, doc);
					
				} else if ("print".equals(command)) {
					printConfigEntryValue(module, name, entryValue);
					
				} else {
					throw new IllegalArgumentException("unknown command: " + command);					
				}
			}

		} catch (Throwable t) {
			t.printStackTrace();
			System.err.println("\nQuitting, no changes written to disk!");
			return;

		}

		// Process changes
		if (!"print".equals(command)) {
			writeChangesToDisk();
		}

	}
	
	private void updateConfigs() throws Exception {
		// update all Chipster configs
		for (String componentDir : getComponentDirsWithConfig()) {
			if (new File(componentDir).exists()) {
				File configFile = new File(componentDir + File.separator + DirectoryLayout.CONF_DIR + File.separator + Configuration.CONFIG_FILENAME);
				updateChipsterConfigFile(configFile, true);
			}
		}
		File wsClientConfigFile = getClientConfigFile();
		if (wsClientConfigFile.exists()) {
			updateChipsterConfigFile(wsClientConfigFile, false);
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
		
		// update CLI client config
		File cliConfigFile = new File(webstartDir + File.separator + DirectoryLayout.WEB_ROOT + File.separator + "chipster-cli.bash");
		if (cliConfigFile.exists()) {
			updateCliClientConfigFile(cliConfigFile);
		}
	}

	private File getClientConfigFile() {
		return new File("webstart" + File.separator + DirectoryLayout.WEB_ROOT + File.separator + Configuration.CONFIG_FILENAME);
	}

	private void updateActivemqConfigFile(File configFile) throws Exception {
		Document doc = openXmlForUpdating("ActiveMQ", configFile, false);
		Element broker = (Element)doc.getDocumentElement().getElementsByTagName("broker").item(0);
		
		Element transportConnectors = (Element)broker.getElementsByTagName("transportConnectors").item(0);		
		Element transportConnector = (Element)transportConnectors.getElementsByTagName("transportConnector").item(0); // edit first in the list (could use attribute name to decide right one..)
		URI uri = new URI(transportConnector.getAttribute("uri"));
		uri.setScheme(configs[BROKER_PROTOCOL_INDEX][VAL_INDEX]);
		uri.setHost(configs[BROKER_PRIVATE_HOST_INDEX][VAL_INDEX]);
		uri.setPort(Integer.parseInt(configs[BROKER_PORT_INDEX][VAL_INDEX]));

		updateElementAttribute(transportConnector, "uri", uri.toString());
		
		writeLater(configFile, doc);
	}
	
	private void updateWsConfigFile(File configFile) throws Exception {
		Document doc = openXmlForUpdating("Web Start", configFile, false);
		Element jnlp = (Element)doc.getDocumentElement();
		updateElementAttribute(jnlp, "codebase", getWebstartUrl());
		Element applicationDesc = (Element)jnlp.getElementsByTagName("application-desc").item(0);
		NodeList arguments = applicationDesc.getElementsByTagName("argument");
		Element lastArgument = (Element)arguments.item(arguments.getLength() - 1);
		String url = "http://" + configs[BROKER_PUBLIC_HOST_INDEX][VAL_INDEX] + ":" + configs[WS_PORT][VAL_INDEX] + "/" + Configuration.CONFIG_FILENAME;
		updateElementValue(lastArgument, "configuration URL (for Web Start)", url);
		writeLater(configFile, doc);
	}
	
	private void updateCliClientConfigFile(File configFile) throws Exception {		
		List<String> conf = openFileForUpdating("Chipster CLI client", configFile, false);
		String hostUrl = "http://" + configs[BROKER_PUBLIC_HOST_INDEX][VAL_INDEX] + ":" + configs[WS_PORT][VAL_INDEX];
		updateStringLine(conf, "host=", "host", "host=\"" + hostUrl + "\"");
		updateStringLine(conf, "conf=", "config", "conf=\"" + Configuration.CONFIG_FILENAME + "\"");			
		writeLater(configFile, conf);
	}

	private void updateActivemqConfigFilePasswords(File configFile) throws Exception {
		Document doc = openXmlForUpdating("ActiveMQ", configFile, false);
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
		Document doc = openXmlForUpdating("Chipster", configFile, false);

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
	
	private void updateChipsterConfigFile(File configFile, boolean isPrivateNetwork) throws Exception {
		Document doc = openXmlForUpdating("Chipster", configFile, false);

		Element messagingModule = XmlUtil.getChildWithAttributeValue(doc.getDocumentElement(), "moduleId", "messaging");
		if (isPrivateNetwork) {
			updateConfigEntryValue(messagingModule, "broker-host", configs[BROKER_PRIVATE_HOST_INDEX][VAL_INDEX]);
		} else {
			updateConfigEntryValue(messagingModule, "broker-host", configs[BROKER_PUBLIC_HOST_INDEX][VAL_INDEX]);			
		}
		updateConfigEntryValue(messagingModule, "broker-protocol", configs[BROKER_PROTOCOL_INDEX][VAL_INDEX]);
		updateConfigEntryValue(messagingModule, "broker-port", configs[BROKER_PORT_INDEX][VAL_INDEX]);

		Element filebrokerModule = XmlUtil.getChildWithAttributeValue(doc.getDocumentElement(), "moduleId", "filebroker");
		if (filebrokerModule != null) {
			updateConfigEntryValue(filebrokerModule, "port", configs[FILEBROKER_PORT_INDEX][VAL_INDEX]);
			updateConfigEntryValue(filebrokerModule, "url", createFilebrokerUrl());
		}

		Element compModule = XmlUtil.getChildWithAttributeValue(doc.getDocumentElement(), "moduleId", "comp");
		if (compModule != null) {
			updateConfigEntryValue(compModule, "max-jobs", configs[MAX_JOBS_INDEX][VAL_INDEX]);
		}
		
		Element webstartModule = XmlUtil.getChildWithAttributeValue(doc.getDocumentElement(), "moduleId", "webstart");
		if (webstartModule != null) {
			updateConfigEntryValue(webstartModule, "port", configs[WS_PORT][VAL_INDEX]);
		}

		Element clientModule = XmlUtil.getChildWithAttributeValue(doc.getDocumentElement(), "moduleId", "client");
		if (clientModule != null) {
			updateConfigEntryValue(clientModule, "manual-root", getWebstartUrl() + "/manual/");
		}
		
		Element managerModule = XmlUtil.getChildWithAttributeValue(doc.getDocumentElement(), "moduleId", "manager");
		if (managerModule != null) {
			updateConfigEntryValue(managerModule, "web-console-port", configs[MANAGER_PORT][VAL_INDEX]);
			updateConfigEntryValue(managerModule, "admin-email", configs[MANAGER_EMAIL][VAL_INDEX]);
		}

		writeLater(configFile, doc);
	}

	private String createFilebrokerUrl() {
		return configs[FILEBROKER_PORT_PROTOCOL][VAL_INDEX] + "://" + configs[BROKER_PUBLIC_HOST_INDEX][VAL_INDEX];
	}

	private void updateConfigEntryValue(Element module, String name, String newValue) {
		Element entry = XmlUtil.getChildWithAttributeValue(module, "entryKey", name);
		Element value = (Element)entry.getElementsByTagName("value").item(0);
		updateElementValue(value, name, newValue);
	}

	private void setConfigEntryValue(Document doc, Element module, String name, String newValue) {
		Element entry = XmlUtil.getChildWithAttributeValue(module, "entryKey", name);
		if (entry == null) {
			entry = doc.createElement("entry");
			module.appendChild(entry);
			entry.setAttribute("entryKey", name);
			Element newValueElement = doc.createElement("value");
			entry.appendChild(newValueElement);
		}
		Element value = (Element)entry.getElementsByTagName("value").item(0);
		updateElementValue(value, name, newValue);
	}

	private void removeConfigEntry(Element module, String name) {
		Element entry = XmlUtil.getChildWithAttributeValue(module, "entryKey", name);
		System.out.println("  removing " + name);
		module.removeChild(entry);
	}

	private void printConfigEntryValue(Element module, String name, String newValue) {
		Element entry = XmlUtil.getChildWithAttributeValue(module, "entryKey", name);
		Element value = (Element)entry.getElementsByTagName("value").item(0);
		System.out.println(value.getTextContent());
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
	
	private void updateStringLine(List<String> lines, String startsWith,
			String logicalName, String newValue) {
		
		for (int i = 0; i < lines.size(); i++) {
			String line = lines.get(i);
			if (line.startsWith(startsWith)) {
				System.out.println("  changing " + logicalName + ": " + line + " -> " + newValue);
				lines.set(i, newValue);
				break;
			}
		}		
	}

	private void writeLater(File configFile, Document doc) throws TransformerException, UnsupportedEncodingException, FileNotFoundException {
		documentsToWrite.put(configFile.getAbsolutePath(), doc);
		System.out.println("");
	}
	

	private void writeLater(File configFile, List<String> conf) {
		filesToWrite.put(configFile.getAbsolutePath(), conf);
		System.out.println("");
	}

	private Document openXmlForUpdating(String name, File configFile, boolean silent) throws SAXException, IOException, ParserConfigurationException {
		if (!silent) {
			System.out.println("Updating " + name + " config in " + configFile.getAbsolutePath());
		}
		Document doc = XmlUtil.parseFile(configFile);
		return doc;
	}
	
	private List<String> openFileForUpdating(String name, File configFile, boolean silent) throws IOException {
		if (!silent) {
			System.out.println("Updating " + name + " config in " + configFile.getAbsolutePath());
		}
		List<String> lines = Files.readAllLines(Paths.get(configFile.toURI()), Charset.defaultCharset());
		return lines;
	}

	private String getWebstartUrl() {
		return "http://" + configs[BROKER_PUBLIC_HOST_INDEX][VAL_INDEX] + ":" + configs[WS_PORT][VAL_INDEX];

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
