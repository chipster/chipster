package fi.csc.microarray.config;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.StringWriter;
import java.net.InetAddress;
import java.util.LinkedList;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

import fi.csc.microarray.util.IOUtils;
import fi.csc.microarray.util.XmlUtil;

/**
 * Simple tool for upgrading existing installation to new version.
 * 
 * @author Aleksi Kallio
 *
 */
public class UpgradeTool {

	private LinkedList<Operation> operations = new LinkedList<Operation>();
	
	private static abstract class Operation {

		private File file;
		private Object parameter;

		public Operation(File file, Object parameter) {
			this.file = file;
			this.parameter = parameter;
		}
		
		public abstract void execute(File file, Object parameter);
		public abstract String describeExecution(File file, Object parameter);
	}
	
	
	public void upgrade(File pathToOld) throws Exception {
		
		// create the transaction
		createOperations(pathToOld);
		
		// tell what we are about to do
		for (Operation operation : operations) {
			System.out.println(operation.describeExecution(operation.file, operation.parameter));
		}
		
		// ask if we should commit
		System.out.println("");
		ConfigTool.verifyChanges(new BufferedReader(new InputStreamReader(System.in)));
		
		// do it
		for (Operation operation : operations) {
			try {
				operation.execute(operation.file, operation.parameter);
				
			} catch (Exception e) {
				e.printStackTrace();
				System.out.println("\nError occurred, continuing anyway...");
			}
		}
		
		// we are done
		System.out.println("\nAll changes successfully written!");
		System.out.println("\nYou need to regenerate server passwords (run genpasswd.sh in the upgraded installation)");

		
	}
	
	public void createOperations(File pathToOld) {

		// transform old installation components to new directory layout
		for (String componentDir : ConfigTool.getComponentDirsWithConfig()) {
			if (new File(pathToOld, componentDir).exists()) {

				// initialise paths
				File compDir = new File(pathToOld, componentDir);
				File workDir12x = new File(compDir, "nami-work-files");
				File logsDir13x = new File(compDir, DirectoryLayout.LOGS_DIR);
				File confDir13x = new File(compDir, DirectoryLayout.CONF_DIR);
				File securityDir13x = new File(compDir, DirectoryLayout.SECURITY_DIR);
				File binDir = new File(compDir, DirectoryLayout.BIN_DIR);

				// transform directory layout
				moveToNewDir(workDir12x, logsDir13x, "messages.log", "jobs.log", "security.log", "status.log");
				delayedMove(new File(workDir12x, "nami.log"), new File(logsDir13x, "chipster.log"));
				moveToNewDir(workDir12x, securityDir13x, "keystore.ks", "users");
				delayedUsersUpgrade(new File(securityDir13x, "users"), new File(workDir12x, "users"));
				confDir13x.mkdirs();
				delayedMove(new File(workDir12x, "jaas.config"), new File(confDir13x, "jaas.config"));
				delayedMove(new File(workDir12x, "nami-config.xml"), new File(confDir13x, "chipster-config.xml"));
				delayedConfigPurge(new File(confDir13x, "chipster-config.xml"), new File(workDir12x, "nami-config.xml"), componentDir);
				delayedMove(new File(compDir, "web-content"), new File(compDir, "web-root"));
				delayedMove(new File(workDir12x, "fileserver"), new File(compDir, "file-root"));
				delayedMove(new File(workDir12x, "analyser-work"), new File(compDir, "jobs-data"));
				delayedDelete(workDir12x);
				
				// keep the later of 32/64 wrapper.conf's
				File[] wrapperConfs = new File[] {
						new File(binDir, "linux-x86-32"  + File.separator + "wrapper.conf"),
						new File(binDir, "linux-x86-64"  + File.separator + "wrapper.conf")
				};	
				
				
				File newWrapperConf = new File(confDir13x, "wrapper.conf");
				if (wrapperConfs[0].lastModified() > wrapperConfs[1].lastModified()) {
					delayedDelete(wrapperConfs[1]);
					delayedMove(wrapperConfs[0], newWrapperConf);					
				} else {					
					delayedDelete(wrapperConfs[0]);
					delayedMove(wrapperConfs[1], newWrapperConf);
				}
				delayedWrapperScriptUpgrade(new File(binDir, "linux-x86-32"  + File.separator + "chipster-" + componentDir));
				delayedWrapperScriptUpgrade(new File(binDir, "linux-x86-64"  + File.separator + "chipster-" + componentDir));				
			}
		}
		
		// replace old stuff in the old installation with new from this one
		File lib = new File(pathToOld, "shared" + File.separator + "lib");
		cleanDir(lib);
		moveAll(new File("shared" + File.separator + "lib"), lib);
		replaceFromNew("client" + File.separator + "bin" + File.separator + "chipster-current.jar", pathToOld);
		delayedMove(new File("webstart" + File.separator + "web-content" + File.separator + "chipster-current.jar"), new File(pathToOld + File.separator + "webstart" + File.separator + "web-root" + File.separator + "chipster-current.jar"));
		replaceFromNew("activemq" + File.separator + "conf" + File.separator + "activemq.xml", pathToOld);
	}

	private void replaceFromNew(String file, File pathToOld) {
		delayedMove(new File(file), new File(pathToOld,  file));
	}

	private void cleanDir(File dir) {
		for (File file : dir.listFiles()) {
			delayedDelete(file);
		}
	}
	
	private void moveAll(File sourceDir, File dest) {
		for (File file : sourceDir.listFiles()) {
			delayedMove(file, new File(dest, file.getName()));
		}
	}

	private void moveToNewDir(File oldDir, File newDir, String... files) {
		
		// create new it does not exist
		newDir.mkdir();
		
		// move listed files to new
		for (String file : files) {
			File from = new File(oldDir, file);
			File to = new File(newDir, file);
			delayedMove(from, to);
		}
	}
	
	private void delayedMove(File from, File to) {

		if (from.exists()) {
			operations.add(new Operation(from, to) {

				public String describeExecution(File file, Object parameter) {
					return "Move " + file + " to " + parameter;
				}

				public void execute(File file, Object parameter) {
					file.renameTo((File)parameter);				
				}			
			});
		}
	}

	private void delayedDelete(File file) {
		if (file.exists()) {
			operations.add(new Operation(file, null) {

				public String describeExecution(File file, Object parameter) {
					return "Delete " + file ;
				}

				public void execute(File file, Object parameter) {
					boolean deleted = file.delete();
					if (!deleted) {
						System.out.println("Warning: could not delete " + file.getAbsolutePath());
					}
				}			
			});
		}
	}
	

	private void delayedUsersUpgrade(File fileAfterUpgrade, File filebeforeUpgrade) {
		if (filebeforeUpgrade.exists()) {
			operations.add(new Operation(fileAfterUpgrade, null) {

				public String describeExecution(File file, Object parameter) {
					return "Updgrade file format in " + file + " (add empty expiration dates if needed)";
				}

				public void execute(File file, Object parameter) {
					BufferedReader in = null;
					OutputStream overwriteOut = null;
					try {
						StringWriter stringWriter = new StringWriter();
						BufferedWriter out = new BufferedWriter(stringWriter);
						in = new BufferedReader(new InputStreamReader(new FileInputStream(file)));
						
						boolean needOverwrite = false;
						
						for (String line = in.readLine(); line != null; line = in.readLine()) {
							if (line.split(":").length == 3) {
								// has exactly two separators, need to add empty for exp. date
								int commentSep = line.lastIndexOf(':');
								String newLine = line.substring(0, commentSep) + ":" + line.substring(commentSep);
								line = newLine;
								needOverwrite = true;
							}
							out.write(line);
						}
						out.close();
						if (needOverwrite) {
							overwriteOut = new FileOutputStream(file);
							IOUtils.copy(new ByteArrayInputStream(stringWriter.getBuffer().toString().getBytes()), overwriteOut);
						}
						
					} catch (IOException e) {
						e.printStackTrace();
						System.out.println("Warning: could not upgrade users file " + file.getAbsolutePath());
						
					} finally {
						IOUtils.closeIfPossible(in);
						IOUtils.closeIfPossible(overwriteOut);
					}
				}			
			});
		}		
	}

	private void delayedConfigPurge(File fileAfterUpgrade, File filebeforeUpgrade, String componentDir) {
		if (filebeforeUpgrade.exists()) {
			operations.add(new Operation(fileAfterUpgrade, componentDir) {

				public String describeExecution(File file, Object parameter) {
					return "Upgrade and purge " + file + " (leave only those modules needed for " + parameter + ")";
				}

				public void execute(File file, Object parameter) {
					BufferedReader in = null;
					OutputStream overwriteOut = null;
					try {
						StringWriter stringWriter = new StringWriter();
						BufferedWriter out = new BufferedWriter(stringWriter);
						in = new BufferedReader(new InputStreamReader(new FileInputStream(file)));
						
						Document originalDocument = XmlUtil.getInstance().parseReader(in);
						IOUtils.closeIfPossible(in);
						Document newDocument = XmlUtil.getInstance().newDocument();
						newDocument.appendChild(newDocument.createElement("configuration"));
						newDocument.getDocumentElement().setAttribute("content-version", "3");

						// pick correct modules
						moveModule("messaging", originalDocument, "messaging", newDocument);
						moveModule("security", originalDocument, "security", newDocument);
						if ("comp".equals(parameter)) {
							moveModule("analyser", originalDocument, "comp", newDocument);
							
						} else if ("frontend".equals(parameter)) {
							moveModule("frontend", originalDocument, "filebroker", newDocument);
							moveModule("filebroker", originalDocument, "filebroker", newDocument);
							
						} else if ("fileserver".equals(parameter)) {
							Element module = (Element)newDocument.getDocumentElement().appendChild(newDocument.createElement("configuration-module"));
							module.setAttribute("moduleId", "filebroker");
							addEntry(module, newDocument, "port", "8080");
							addEntry(module, newDocument, "url", "http://" + InetAddress.getLocalHost().getHostName());
							
							
						} else if ("webstart".equals(parameter)) {
							Element module = (Element)newDocument.getDocumentElement().appendChild(newDocument.createElement("configuration-module"));
							module.setAttribute("moduleId", "webstart");
							addEntry(module, newDocument, "port", "8081");
							
						} else if ("manager".equals(parameter)) {
							Element module = (Element)newDocument.getDocumentElement().appendChild(newDocument.createElement("configuration-module"));
							module.setAttribute("moduleId", "manager");
							addEntry(module, newDocument, "port", "8082");

						} else {
							moveModuleIfExists((String)parameter, originalDocument, (String)parameter, newDocument);
						}
						
						// remove obsolete entries
						NodeList entries = newDocument.getElementsByTagName("entry");
						for (int i = 0; i < entries.getLength(); i++) {
							Element entry = (Element)entries.item(i);
							String oldName = entry.getAttribute("entryKey");
							
							if ("filebroker_urls".equals(oldName)) {
								entry.getParentNode().removeChild(entry);
								break;
							}
						}
						
						// fix entry names
						entries = newDocument.getElementsByTagName("entry");
						for (int i = 0; i < entries.getLength(); i++) {
							Element entry = (Element)entries.item(i);
							String oldName = entry.getAttribute("entryKey");

							// rename
							oldName = oldName.replace('_', '-');
							String newName = "";
							for (int c = 0; c < oldName.length(); c++) {
								char character = oldName.charAt(c);
								if (Character.isUpperCase(character)) {
									newName += '-';
									newName += Character.toLowerCase(character);
								} else {
									newName += character;
								}
							}
							
							if (newName.startsWith("-")) {
								newName = newName.substring(1);
							}
							
							entry.setAttribute("entryKey", newName);
							
						}
												
						// overwrite previous
						XmlUtil.getInstance().printXml(newDocument, out);
						out.close();

						overwriteOut = new FileOutputStream(file);
						IOUtils.copy(new ByteArrayInputStream(stringWriter.getBuffer().toString().getBytes()), overwriteOut);
						
					} catch (Exception e) {
						e.printStackTrace();
						System.out.println("Warning: could not upgrade configuration file " + file.getAbsolutePath());
						
					} finally {						
						IOUtils.closeIfPossible(overwriteOut);
					}
				}

				private void addEntry(Element module, Document document, String name, String value) {
					Element entry = (Element)module.appendChild(document.createElement("entry"));
					entry.setAttribute("entryKey", name);
					Element valueElement = (Element)entry.appendChild(document.createElement("value"));
					valueElement.setTextContent(value);
				}

				private void moveModuleIfExists(String originalModule, Document originalDocument, String newModule, Document newDocument) {
					try {
						moveModule(originalModule, originalDocument, newModule, newDocument);
					} catch (IllegalArgumentException e) {
						// module did not exist, ignore the exception 
					}
				}
				private void moveModule(String originalModule, Document originalDocument, String newModule, Document newDocument) {
					NodeList modules = originalDocument.getElementsByTagName("configuration-module");
					for (int i = 0; i < modules.getLength(); i++) {
						Element module = (Element)modules.item(i);
						if (originalModule.equals(module.getAttribute("moduleId"))) {
							Element importedElement = (Element)newDocument.importNode(module, true);
							newDocument.getDocumentElement().appendChild(importedElement);
							importedElement.setAttribute("moduleId", newModule);
							return;
						}
					}
					throw new IllegalArgumentException("module " + originalModule + " was not found from config");
				}			
			});
		}		
	}
	
	private void delayedWrapperScriptUpgrade(File file) {
		if (file.exists()) {
			operations.add(new Operation(file, null) {

				public String describeExecution(File file, Object parameter) {
					return "Updgrade " + file + " (change wrapper.conf path)";
				}

				public void execute(File file, Object parameter) {
					OutputStream overwriteOut = null;
					try {
						String wrapperScript = fileToString(file);
						
						wrapperScript = wrapperScript.replace("./wrapper.conf", "../../conf/wrapper.conf");
						
						overwriteOut = new FileOutputStream(file);						
						IOUtils.copy(new ByteArrayInputStream(wrapperScript.getBytes()), overwriteOut);
						
					} catch (IOException e) {
						e.printStackTrace();
						System.out.println("Warning: could not upgrade wrapper script file " + file.getAbsolutePath());
						
					} finally {
						IOUtils.closeIfPossible(overwriteOut);
					}
				}
				
				public String fileToString(File file) throws IOException  {
					StringBuffer buffer = new StringBuffer();
					BufferedReader inputReader = new BufferedReader(new InputStreamReader(new FileInputStream(file)));
					String line;
					for (line = inputReader.readLine(); line != null; line = inputReader.readLine()) {
						buffer.append(line + "\n");
					}
					
					return buffer.toString();
				}

			});
		}		
	}

}


