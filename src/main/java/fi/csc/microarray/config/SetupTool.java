package fi.csc.microarray.config;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.net.MalformedURLException;
import java.net.URL;
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
 * Simple tool for setting up external environment for Chipster comp service. This
 * tools install OS packages, applications (download & compile) and various types of 
 * R packages.
 * 
 * @author Aleksi Kallio
 *
 */
public class SetupTool {

	private static interface Installer {
		public boolean install(Element item, Element dependsOnBundle) throws Exception;
	}
	
	//private final File ENVIRONMENT_XML = new File("comp/conf/environment.xml");
	private final File ENVIRONMENT_XML = new File("/home/akallio/eclipse-workspace/chipster/src/main/applications/wrapper/comp/conf/environment.xml");
	
	public static void main(String[] args) throws Exception {
		new SetupTool().setup("");
	}
	
	public void setup(String command) throws Exception {
		
		// initialise the setup... or setup the initialisation?
		System.out.println("Chipster environment installer processing " + ENVIRONMENT_XML.getAbsolutePath());
		HashMap<String, Installer> installers = new HashMap<String, Installer>();
		installers.put("application", new ApplicationInstaller());
		installers.put("r-package", new RPackageInstaller());
		installers.put("bioconductor", new BioconductorInstaller());
		installers.put("os-package", new OSPackageInstaller());
		
		Document xml = XmlUtil.parseFile(ENVIRONMENT_XML);
		
		// iterate over installation bundles and gather work items
		NodeList tasks = xml.getDocumentElement().getElementsByTagName("bundle");
		Queue<Element> workQueue = new LinkedList<Element>();		
		
		for (int i = 0; i < tasks.getLength(); i++) {
			Element bundle = (Element)tasks.item(i); 
			NodeList children = bundle.getChildNodes();
			for (int j = 0; j < children.getLength(); j++) {
				if (children.item(j) instanceof Element ) {
					Element element = (Element)children.item(j);
					if (isInstalled(element)) {
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
				String depends = ((Element)e.getParentNode()).getAttribute("depends");
				if (!"".equals(depends)) {
					requiredBundle = XmlUtil.getChildWithAttributeValue(xml.getDocumentElement(), "name", depends);
					
					if (!isInstalled(requiredBundle)) {
						// can not install at least yet, because dependences are not installed
						continue;
					}
				}
				
				// found something to work with
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

				Installer installer = installers.get(element.getAttribute("type"));
				if (installer == null) {
					System.out.println("unrecognised item type: " + element.getAttribute("type"));
					failedItems.add(element);
					continue;
				}

				boolean ok = installer.install(element, requiredBundle);

				// mark this done if it went okay
				if (ok) {
					System.out.println("  Installation succesfull.");
					element.setAttribute("installed", "yes");
					successfullItems.add(element);

				} else {
					System.out.println("  Installation FAILED.");
					failedItems.add(element);
				}				
			}						
		}		
		
		// finally check how well the process went 
		System.out.println("\n\n********************\nINSTALLATION SUMMARY\n********************\n");
		       
		if (!successfullItems.isEmpty()) {
			System.out.println("Installation for following items was successfull");
			for (Element e : successfullItems) {
				System.out.println("  " + descibeItem(e));
			}			
		}				
		boolean requiredNotInstalled = false;
		if (!failedItems.isEmpty()) {
			System.out.println("Installation for following items failed");
			for (Element e : failedItems) {
				System.out.println("  " + descibeItem(e));
				if ("true".equals(((Element)e.getParentNode()).getAttribute("required"))) {
					requiredNotInstalled = true;
				}
			}			
		}		
		if (!workQueue.isEmpty()) {
			System.out.println("\nFollowing items could not be installed because of unmet dependencies");
			for (Element e : workQueue) {
				System.out.println("  " + descibeItem(e));
				if ("true".equals(((Element)e.getParentNode()).getAttribute("required"))) {
					requiredNotInstalled = true;
				}
			}
		}

		// print out last summary line
		System.out.print("\nResult: ");
		if (requiredNotInstalled) {
			System.out.println("Required bundles were not installed. Environment is not usable.");
		} else if (!workQueue.isEmpty() || !failedItems.isEmpty() ) {
			System.out.println("Some bundles were not installed, but all required bundles were. Environment is usable but not complete. ");
		} else {
			System.out.println("All bundles were installed. Environment is completely installed.");
		}
	}
	
	private boolean isInstalled(Element element) {
		if ("item".equals(element.getLocalName())) {
			return "true".equals(element.getAttribute("installed"));
			
		} else if ("bundle".equals(element.getLocalName())) {
			NodeList items = element.getElementsByTagName("item");
			for (int i = 0; i < items.getLength(); i++) {
				if (!isInstalled((Element)items.item(i))) {
					// one item inside bundle not installed => bundle not installed
					return false; 
				}
			}
			
			// everything was installed inside the bundle
			return true;
			
		} else {
			throw new RuntimeException("illegal XML element: " + element.getLocalName());
		}
		
	}
	
	private String descibeItem(Element item) {
		return item.getLocalName() + " " + item.getAttribute("name") + " (in bundle " + ((Element)item.getParentNode()).getAttribute("name") + ")";	
	}
	
	private static class ApplicationInstaller implements Installer {

		public boolean install(Element item, Element dependsOnBundle) {
			return true; // false
		}
		
	}
	
	private static class RPackageInstaller implements Installer {

		public boolean install(Element item, Element dependsOnBundle) throws Exception {
			try {
				Element rBundle = XmlUtil.getChildWithAttributeValue(dependsOnBundle, "type", "application");
				String rExecutable = ((Element)rBundle.getElementsByTagName("dir").item(0)).getAttribute("path") + "/bin/R";
				rExecutable = rExecutable.replaceAll("//", "/");

				// see how this package is available
				if (item.getElementsByTagName("url").getLength() > 0) {
					// available as .tar.gz via URL, download it
					URL url = new URL(item.getElementsByTagName("url").item(0).getTextContent().trim());
					return downloadAndInstall(url, rExecutable);
					
				} else if (item.getElementsByTagName("repository").getLength() > 0) {
					String repository = item.getElementsByTagName("repository").item(0).getTextContent().trim();
					String packageName = item.getElementsByTagName("package").item(0).getTextContent().trim();
					return runRCommand(rExecutable, "install.packages(c(\"" + packageName + "\"), repos=\"" + repository + "\", dependencies = T)");

				} else if (item.getElementsByTagName("default-bioconductor-packages").getLength() > 0) {
					return runRCommands(rExecutable, new String[] {
							"source(\"http://www.bioconductor.org/biocLite.R\")",
							"biocLite()"
					});

				} else if (item.getElementsByTagName("bioconductor-package").getLength() > 0) {
					String packageName = item.getElementsByTagName("bioconductor-package").item(0).getTextContent().trim();
					return runRCommands(rExecutable, new String[] {
							"source(\"http://www.bioconductor.org/biocLite.R\")",
							"biocLite(c(\"" + packageName + "\"))"
					});

				} else if (item.getElementsByTagName("from-web-page").getLength() > 0) {
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

		private boolean downloadAndInstall(URL url, String rExecutable) throws MalformedURLException, JMSException, IOException, FileNotFoundException, InterruptedException {
			System.out.println("  Download and install package from " + url + " to " + rExecutable);
			File file = downloadFile(url);				

			// install the package
			boolean ok = runRCommand(rExecutable, "install.packages(pkgs=\"" + file.getAbsolutePath() + "\", repos=NULL)");
			
			// clean up
			file.delete();
			return ok;
		}

		private boolean runRCommand(String rExecutable, String command) throws IOException, InterruptedException {
			return runRCommands(rExecutable, new String[] { command} );
		}
		private boolean runRCommands(String rExecutable, String[] commands) throws IOException, InterruptedException {
			Process process = Runtime.getRuntime().exec((rExecutable + " --vanilla"), null, new File(System.getProperty("user.dir")));
			ByteArrayOutputStream bufferOut = new ByteArrayOutputStream();
			ByteArrayOutputStream bufferErr = new ByteArrayOutputStream();
			startBackgroundEchoThread(process.getInputStream(), new OutputStream[] { System.out, bufferOut } );
			startBackgroundEchoThread(process.getErrorStream(), new OutputStream[] { System.err, bufferErr });			
			PrintWriter commandWriter = new PrintWriter(process.getOutputStream());
			for (String command : commands) {
				commandWriter.println(command);
			}
			commandWriter.println("quit(save=\"no\")");
			commandWriter.flush();
			process.waitFor();
			
			return bufferOut.toString().contains("DONE");
		}

		private void startBackgroundEchoThread(final InputStream inputStream, final OutputStream[] outs) {
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

		private File downloadFile(URL url) throws JMSException, IOException, MalformedURLException, FileNotFoundException {
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
		
		private String downloadString(URL url) throws JMSException, IOException, MalformedURLException, FileNotFoundException {
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
	
	private static class BioconductorInstaller implements Installer {

		public boolean install(Element item, Element dependsOnBundle) {
			return false;
		}
		
		
	}
	private static class OSPackageInstaller implements Installer {

		public boolean install(Element item, Element dependsOnBundle) {
			return true; // false
		}
		
	}
	
}


