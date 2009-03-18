package fi.csc.microarray.config;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.util.LinkedList;

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
				delayedMove(new File(workDir12x, "jaas.config"), new File(confDir13x, "jaas.config"));
				delayedMove(new File(workDir12x, "nami-config.xml"), new File(confDir13x, "chipster-config.xml"));
				delayedMove(new File(compDir, "web-content"), new File(compDir, "web-root"));
				delayedMove(new File(workDir12x, "fileserver"), new File(compDir, "file-root"));
				delayedMove(new File(workDir12x, "analyser-work"), new File(compDir, "jobs-data"));
				delayedDelete(workDir12x);
				
				// keep the later of 32/64 wrapper.conf's
				File[] wrapperConfs = new File[] {
						new File(binDir, "linux-x86-32"  + File.separator + "wrapper.conf"),
						new File(binDir, "linux-x86-64"  + File.separator + "wrapper.conf")
				};				
				if (wrapperConfs[0].lastModified() > wrapperConfs[1].lastModified()) {
					delayedMove(wrapperConfs[0], confDir13x);
					delayedDelete(wrapperConfs[1]);
				} else {					
					delayedDelete(wrapperConfs[0]);
					delayedMove(wrapperConfs[1], confDir13x);
				}
				
				// update file content to new version
				delayedConfigPurge(new File(confDir13x, "chipster-config.xml"), componentDir);
				delayedUsersUpgrade(new File(securityDir13x, "users"));
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
						System.out.println("Warning: could not delete " + file);
					}
				}			
			});
		}
	}
	

	private void delayedUsersUpgrade(File file) {
		if (file.exists()) {
			operations.add(new Operation(file, null) {

				public String describeExecution(File file, Object parameter) {
					return "Updgrade file format in " + file + " (add empty expiration dates if needed)";
				}

				public void execute(File file, Object parameter) {
					// TODO				
				}			
			});
		}		
	}

	private void delayedConfigPurge(File file, String componentDir) {
		if (file.exists()) {
			operations.add(new Operation(file, componentDir) {

				public String describeExecution(File file, Object parameter) {
					return "Purge unneeded configuration modules from " + file + " (leave only those needed for " + parameter + ")";
				}

				public void execute(File file, Object parameter) {
					// TODO				
				}			
			});
		}		
	}
}


