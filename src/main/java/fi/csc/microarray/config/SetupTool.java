package fi.csc.microarray.config;

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.CharArrayWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.net.MalformedURLException;
import java.net.URL;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Queue;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.jms.JMSException;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

import fi.csc.microarray.util.IOUtils;
import fi.csc.microarray.util.UrlTransferUtil;
import fi.csc.microarray.util.XmlUtil;

/**
 * Simple tool for setting up external environment for Chipster comp service. This tools install OS packages, applications (download &
 * compile) and various types of R packages. In version 1.4.0 only R packages are supported.
 * 
 * @author Aleksi Kallio
 * 
 */

public class SetupTool {

	private static interface Installer {
		public boolean install(Element item, Element dependsOnBundle, PrintWriter infoWriter) throws Exception;
	}

	private final File ENVIRONMENT_XML = new File("comp/conf/environment.xml");
	
	public static int main(String[] args) throws Exception {
		return new SetupTool().setup();
	}

	public int setup() throws Exception {

		// start
		BufferedReader in = new BufferedReader(new InputStreamReader(System.in)); 
		System.out.println("Chipster environment installer processing " + ENVIRONMENT_XML.getAbsolutePath());
		System.out.println("Please press any key to continue or Ctrl+C to abort...");
		in.readLine();
		
		// initialise the setup... or setup the initialisation?
		HashMap<String, Installer> installers = new HashMap<String, Installer>();
		installers.put("application", new DummyInstaller());
		installers.put("file", new DummyInstaller());
		installers.put("r-package", new RPackageInstaller());
		installers.put("os-package", new OSPackageInstaller());

		// initialise logging
		String date = new SimpleDateFormat("yyyyMMdd").format(new Date());
		String time = new SimpleDateFormat("kkmm").format(new Date());
		File logFile = createLogFile("install", date, time);
		File infoFile = createLogFile("extra.info", date, time);
		PrintWriter logOut = null;
		PrintWriter infoOut = null;
		FileOutputStream xmlOut = null;
		
		try {
			logOut = new PrintWriter(new FileOutputStream(logFile));
			infoOut = new PrintWriter(new FileOutputStream(infoFile));

			// read environment description
			Document xml = XmlUtil.parseFile(ENVIRONMENT_XML);

			// iterate over installation bundles and gather work items
			NodeList tasks = xml.getDocumentElement().getElementsByTagName("bundle");
			Queue<Element> workQueue = new LinkedList<Element>();

			for (int i = 0; i < tasks.getLength(); i++) {
				Element bundle = (Element) tasks.item(i);
				NodeList children = bundle.getChildNodes();
				for (int j = 0; j < children.getLength(); j++) {
					if (children.item(j) instanceof Element) {
						Element element = (Element) children.item(j);
						if (isInstalled(element)) {
							logOut.println("Already installed: " + descibeItem(element));
							continue; // already installed, skip
						}
						workQueue.add(element);
					}
				}

			}

			// process work items
			Queue<Element> failedItems = new LinkedList<Element>();
			LinkedList<Element> successfullItems = new LinkedList<Element>();
			while (!workQueue.isEmpty()) {

				// search something to work with
				Element element = null;
				Element requiredBundle = null;
				for (Element e : workQueue) {
					String depends = ((Element) e.getParentNode()).getAttribute("depends");
					if (!"".equals(depends)) {
						requiredBundle = XmlUtil.getChildWithAttributeValue(xml.getDocumentElement(), "name", depends);

						if (!isInstalled(requiredBundle)) {
							// can not install at least yet, because dependences are not installed
							continue;
						}
					}

					// found something to work on
					element = e;
					workQueue.remove(e);
					break;
				}

				// if something was found, process it
				if (element == null) {
					// nothing was found
					break;

				} else {
					System.out.println("\nInstalling " + descibeItem(element) + "...");

					// find correct installer for this item
					Installer installer = installers.get(element.getAttribute("type"));
					if (installer == null) {
						System.out.println("unrecognised item type: " + element.getAttribute("type"));
						failedItems.add(element);
						continue;
					}

					// do the actual work
					boolean ok = installer.install(element, requiredBundle, infoOut);

					// mark this done if it went okay
					if (ok) {
						System.out.println("  Installation succesfull.");
						logOut.println("Installed successfully: " + descibeItem(element));
						element.setAttribute("installed", "true");
						successfullItems.add(element);

					} else {
						System.out.println("  Installation FAILED.");
						logOut.println("Installation FAILED: " + descibeItem(element));
						failedItems.add(element);
					}
				}
			}

			// finally check how well the process went
			CharArrayWriter caout = new CharArrayWriter();
			PrintWriter out = new PrintWriter(caout);
			out.println("\n\n********************\nINSTALLATION SUMMARY\n********************\n");

			if (!successfullItems.isEmpty()) {
				out.println("Installation for following items was successfull");
				for (Element e : successfullItems) {
					out.println("  " + descibeItem(e));
				}
			}
			boolean requiredNotInstalled = false;
			if (!failedItems.isEmpty()) {
				out.println("Installation for following items failed");
				for (Element e : failedItems) {
					out.println("  " + descibeItem(e));
					if ("true".equals(((Element) e.getParentNode()).getAttribute("required"))) {
						requiredNotInstalled = true;
					}
				}
			}
			if (!workQueue.isEmpty()) {
				out.println("\nFollowing items could not be installed because of unmet dependencies");
				for (Element e : workQueue) {
					out.println("  " + descibeItem(e));
					if ("true".equals(((Element) e.getParentNode()).getAttribute("required"))) {
						requiredNotInstalled = true;
					}
				}
			}

			// decide the final result
			int returnCode;
			out.print("\nResult: ");
			if (requiredNotInstalled) {
				out.println("Required bundles were not installed. Environment is not usable.");
				returnCode = 2; // full fail
			} else if (!workQueue.isEmpty() || !failedItems.isEmpty()) {
				out.println("Some bundles were not installed, but all required bundles were. Environment is usable but not complete. ");
				returnCode = 1; // partial fail
			} else {
				out.println("All bundles were installed. Environment is completely installed.");
				returnCode = 0; // ok
			}

			// write buffer out
			out.close();
			System.out.println(caout.toString());
			System.out.println("Written installation log to " + logFile);
			System.out.println("Written additional instructions for required manual installation steps (if any) to " + infoFile);

			// write changed installation states
			xmlOut = new FileOutputStream(ENVIRONMENT_XML);
			XmlUtil.printXml(xml, xmlOut);
			
			return returnCode;
			
		} finally {
			IOUtils.closeIfPossible(infoOut);
			IOUtils.closeIfPossible(logOut);
			IOUtils.closeIfPossible(xmlOut);
		}
	}

	private File createLogFile(String name, String date, String time) {
		File file = new File(name + "." + date + "." + time + ".log");
		int i = 2;
		while (file.exists()) {
			file = new File(name + "." + date + "." + time + "." + i + ".log");
			i++;
		}
		return file;
	}

	private boolean isInstalled(Element element) {
		if ("item".equals(element.getLocalName())) {
			return "true".equals(element.getAttribute("installed"));

		} else if ("bundle".equals(element.getLocalName())) {
			NodeList items = element.getElementsByTagName("item");
			for (int i = 0; i < items.getLength(); i++) {
				if (!isInstalled((Element) items.item(i))) {
					// one item inside bundle not installed => bundle not installed
					return false;
				}
			}

			// all items inside the bundle were installed => bundle is installed
			return true;

		} else {
			throw new RuntimeException("illegal XML element: " + element.getLocalName());
		}

	}

	private String descibeItem(Element item) {
		return item.getLocalName() + " " + item.getAttribute("name") + " (in bundle " + ((Element) item.getParentNode()).getAttribute("name") + ")";
	}

	private static class DummyInstaller implements Installer {

		public boolean install(Element item, Element dependsOnBundle, PrintWriter infoWriter) {

			// You should install and compile applications here. Now we just print out instructions			
			for (Element command : XmlUtil.getChildElements(item)) {
				if ("dir".equals(command.getNodeName())) {
					infoWriter.println("In " + command.getTextContent() + ": ");
					
				} else if ("comment".equals(command.getNodeName())) {
					String comment = command.getTextContent().trim();
					for (String line : comment.split("\n")) {
						infoWriter.println("  " + line.trim());
					}
					
					
				} else {
					throw new RuntimeException("unknown element in application item: " + command.getNodeName());
				}
			}
			infoWriter.println("");
			
			return true;
		}

	}

	private static class RPackageInstaller implements Installer {

		public boolean install(Element item, Element dependsOnBundle, PrintWriter infoWriter) throws Exception {
			try {
				String rExecutable = fetchRExecutable(dependsOnBundle);

				// see how this package is available
				if (item.getElementsByTagName("url").getLength() > 0) {
					// available as .tar.gz via URL, download it
					URL url = new URL(item.getElementsByTagName("url").item(0).getTextContent().trim());
					return downloadAndInstall(url, rExecutable);

				} else if (item.getElementsByTagName("repository").getLength() > 0) {
					// available in CRAN repository
					String repository = item.getElementsByTagName("repository").item(0).getTextContent().trim();
					String packageName = item.getElementsByTagName("package").item(0).getTextContent().trim();
					if (isRPackageInstalled(rExecutable, packageName)) {
						return true; // can be skipped
					}
					return runRInstallCommand(rExecutable, "install.packages(c(\"" + packageName + "\"), repos=\"" + repository + "\")");

				} else if (item.getElementsByTagName("default-bioconductor-packages").getLength() > 0) {
					// default BioC Lite packages available via Bioconductor installation mechanism
					String mirror = item.getElementsByTagName("mirror").item(0).getTextContent().trim();
					return runRInstallCommands(rExecutable, generateBiocCommands(new String[] { 
							"biocLite()" 
					}, mirror));

				} else if (item.getElementsByTagName("bioconductor-package").getLength() > 0) {
					// available via Bioconductor installation mechanism
					String packageName = item.getElementsByTagName("bioconductor-package").item(0).getTextContent().trim();
					String mirror = item.getElementsByTagName("mirror").item(0).getTextContent().trim();
					if (isRPackageInstalled(rExecutable, packageName)) {
						return true; // can be skipped
					}
					return runRInstallCommands(rExecutable, generateBiocCommands(new String[] { 
							"biocLite(c(\"" + packageName + "\"))" 
					}, mirror));

				} else if (item.getElementsByTagName("bioconductor-repository").getLength() > 0) {
					// Bioconductor annotation repository
					String repositoryName = item.getElementsByTagName("bioconductor-repository").item(0).getTextContent().trim();
					String mirror = item.getElementsByTagName("mirror").item(0).getTextContent().trim();
					return runRInstallCommands(rExecutable, generateBiocCommands(new String[] { 
							"setRepositories(ind=c(" + repositoryName + "))",
							"install.packages(available.packages())"
					}, mirror));

				} else if (item.getElementsByTagName("from-web-page").getLength() > 0) {
					// available as .tar.gz via URL's listed on web page => web crawl them
					URL url = new URL(item.getElementsByTagName("from-web-page").item(0).getTextContent().trim());
					String pageContents = downloadString(url);

					// match URL's inside <a href="...">
					Matcher matcher = Pattern.compile("href=\"([^\"]*)\"").matcher(pageContents);

					// pick links to R source packages
					while (matcher.find()) {
						String link = matcher.group(1);
						if (link.endsWith(".tar.gz")) {

							if (!link.startsWith("http")) {
								// relative URL
								link = url.toString().substring(0, url.toString().lastIndexOf('/') + 1) + link;
							}

							// install this package
							boolean ok = downloadAndInstall(new URL(link), rExecutable);
							if (!ok) {
								return false;
							}
						}
					}

					return true;

				} else {
					throw new RuntimeException("don't know how to run " + item.getAttribute("name"));
				}

			} catch (Exception e) {
				e.printStackTrace();
				return false;
			}
		}

		private String[] generateBiocCommands(String[] commands, String mirror) {
			String[] initBioc = new String[] {
					"source(\"http://www.bioconductor.org/biocLite.R\")",
					"options(\"BioC_mirror\" = c(\"Mirror\"=\"" + mirror + "\"))"
			};

			String[] biocCommands = new String[initBioc.length + commands.length];
			System.arraycopy(initBioc, 0, biocCommands, 0, initBioc.length);
			System.arraycopy(commands, 0, biocCommands, initBioc.length, commands.length);
			
			return biocCommands;
		}

		private String fetchRExecutable(Element dependsOnBundle) {
			Element rBundle = XmlUtil.getChildWithAttributeValue(dependsOnBundle, "type", "application");
			String rExecutable = ((Element) rBundle.getElementsByTagName("dir").item(0)).getTextContent() + "/bin/R";
			rExecutable = rExecutable.replaceAll("//", "/");
			return rExecutable;
		}

		private boolean downloadAndInstall(URL url, String rExecutable) throws MalformedURLException, JMSException, IOException, FileNotFoundException, InterruptedException {
			System.out.println("  Download and install package from " + url + " to " + rExecutable);
			File file = downloadFile(url);

			// install the package
			boolean ok = runRInstallCommand(rExecutable, "install.packages(pkgs=\"" + file.getAbsolutePath() + "\", repos=NULL)");

			// clean up
			file.delete();
			return ok;
		}
	}

	private static class OSPackageInstaller implements Installer {

		
		public boolean install(Element item, Element dependsOnBundle, PrintWriter infoWriter) {
			// Should install OS packages, but now we just print out instructions

			String name = XmlUtil.getChildElement(item, "package-name").getTextContent();
			Element altName = XmlUtil.getChildElement(item, "alternative-package-name");
			if (altName == null) {
				infoWriter.println("Install manually OS package " + name);
			} else {
				infoWriter.println("Install manually OS package " + name + " (or alternatively " + altName.getTextContent() + ")");
			}
			infoWriter.println("");
			
			return true;
		}

	}

	
	private static boolean isRPackageInstalled(String rExecutable, String packageName) throws IOException, InterruptedException {
		// in the R code we do some tricks because the source code itself is echoed to output 
		String buffer = runRCommands(rExecutable, new String[] { "if(require(\"" + packageName + "\")) { print(gsub(\"_remove this_\", \"\", paste(\"IS ALREADY INSTALLED:_remove this_\", \"" + packageName + "\"))) }" });
		return buffer.contains("IS ALREADY INSTALLED: " + packageName);
	}
	
	private static boolean runRInstallCommand(String rExecutable, String command) throws IOException, InterruptedException {
		return runRInstallCommands(rExecutable, new String[] { command });
	}

	private static boolean runRInstallCommands(String rExecutable, String[] commands) throws IOException, InterruptedException {
		String buffer =  runRCommands(rExecutable, commands);

		// check that something was done and nothing failed (in case multiple packages installed due to dependencies)
		return buffer.contains("DONE") && !buffer.contains("FAILED") && !buffer.contains("ERROR");
	}
	
	private static String runRCommands(String rExecutable, String[] commands) throws IOException, InterruptedException {
		Process process = Runtime.getRuntime().exec((rExecutable + " --vanilla"), null, new File(System.getProperty("user.dir")));
		ByteArrayOutputStream buffer = new ByteArrayOutputStream();
		startBackgroundEchoThread(process.getInputStream(), new OutputStream[] { System.out, buffer });
		startBackgroundEchoThread(process.getErrorStream(), new OutputStream[] { System.err, buffer });
		PrintWriter commandWriter = new PrintWriter(process.getOutputStream());
		for (String command : commands) {
			commandWriter.println(command);
		}
		commandWriter.println("quit(save=\"no\")");
		commandWriter.flush();
		process.waitFor();

		return buffer.toString(); 
	}

	private static void startBackgroundEchoThread(final InputStream inputStream, final OutputStream[] outs) {
		new Thread(new Runnable() {
			public void run() {
				try {
					for (int i = inputStream.read(); i != -1; i = inputStream.read()) {
						for (OutputStream out : outs) {
							out.write(i);
						}
					}
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}).start();
	}

	private static File downloadFile(URL url) throws JMSException, IOException, MalformedURLException, FileNotFoundException {
		File file = new File(new File(url.getFile()).getName());
		InputStream downloadStream = null;
		FileOutputStream out = null;
		file.deleteOnExit();
		try {
			downloadStream = UrlTransferUtil.downloadStream(url);
			out = new FileOutputStream(file);
			IOUtils.copy(downloadStream, out);

		} finally {
			IOUtils.closeIfPossible(downloadStream);
			IOUtils.closeIfPossible(out);
		}
		return file;
	}

	private static String downloadString(URL url) throws JMSException, IOException, MalformedURLException, FileNotFoundException {
		InputStream downloadStream = null;
		ByteArrayOutputStream out = new ByteArrayOutputStream();

		try {
			downloadStream = UrlTransferUtil.downloadStream(url);
			IOUtils.copy(downloadStream, out);

		} finally {
			IOUtils.closeIfPossible(downloadStream);
			IOUtils.closeIfPossible(out);
		}
		return out.toString();
	}

}
