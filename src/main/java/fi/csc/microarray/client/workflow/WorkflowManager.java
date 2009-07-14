package fi.csc.microarray.client.workflow;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.LinkedList;
import java.util.List;

import bsh.EvalError;
import bsh.Interpreter;
import fi.csc.microarray.client.AtEndListener;
import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.dialog.DialogInfo.Severity;
import fi.csc.microarray.client.workflow.api.WfApplication;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.util.Exceptions;
import fi.csc.microarray.util.GeneralFileFilter;


/**
 * 
 * @author Aleksi Kallio
 *
 */
public class WorkflowManager {
	
	
	public static final String WORKFLOW_VERSION = "BSH/2";

	public static void checkVersionHeaderLine(String line) throws IllegalArgumentException {
		if (!line.contains(WORKFLOW_VERSION + " ")) {
			throw new IllegalArgumentException("Script version not supported. Supported version is " + WORKFLOW_VERSION + ", but script begins with " + line);
		}
	}

	public static final String SCRIPT_EXTENSION = "bsh";
	String[] extensions = { WorkflowManager.SCRIPT_EXTENSION };
	public static final GeneralFileFilter FILE_FILTER = 
		new GeneralFileFilter("Workflow in BeanShell format", new String[]{ WorkflowManager.SCRIPT_EXTENSION });
	
	public File scriptDirectory;

	private ClientApplication application;

	public WorkflowManager(ClientApplication application) throws IOException {
		this.application = application;
		this.scriptDirectory = DirectoryLayout.getInstance().getUserDataDir();
	}

	public List<File> getWorkflows() {
		LinkedList<File> workflows = new LinkedList<File>();
		File[] scripts = scriptDirectory.listFiles(new FilenameFilter() {
			public boolean accept(File dir, String name) {
				return name.endsWith(SCRIPT_EXTENSION);
			}				
		});
		for (File script : scripts) {
			workflows.add(script);
		}
		return workflows;
	}

	public void runScript(final File file, final AtEndListener listener) {
		try {
			runScript(file.toURL(), listener);			
		} catch (MalformedURLException e) {			
			throw new IllegalArgumentException(e);
		}
	}
	
	public void runScript(final URL workflowUrl, final AtEndListener listener) {

		Runnable runnable = new Runnable() {
			public void run() {
				BufferedReader in = null;
				try {
					
					in = new BufferedReader(new InputStreamReader(workflowUrl.openConnection().getInputStream()));
					String line = in.readLine();
					checkVersionHeaderLine(line);
					in.close();
					Interpreter i = initialiseBshEnvironment();
					i.eval(new InputStreamReader(workflowUrl.openConnection().getInputStream()));
					
				} catch (Throwable e) {
					e.printStackTrace();
					application.showDialog("Running workflow failed", "Running workflow " + workflowUrl.getFile() + " failed. Usually this is because the workflow contained operations that do not work with the data currently in use. For more information please see details below.", Exceptions.getStackTrace(e), Severity.WARNING, true);
					
				} finally {
					if (listener != null) {
						listener.atEnd();
					}
					if (in != null) {
						try {
							in.close(); // might be closed twice but that is legal
						} catch (IOException e) {} 
					}
				}
			}
		};
		
		new Thread(runnable).start();
	}

	public Interpreter initialiseBshEnvironment() {
		Interpreter i = new Interpreter();
		try {
			i.set("app", new WfApplication(Session.getSession().getApplication()));
		} catch (EvalError ee) {
			throw new RuntimeException("BeanShell console failed to open: " + ee.getMessage());
		}
		return i;	
	}

	public void saveSelectedWorkflow(File selectedFile) throws IOException {
		WorkflowWriter writer = new WorkflowWriter();
		StringBuffer script = writer.writeWorkflow(application.getSelectionManager().getSelectedDataBean());
		saveScript(selectedFile, script);
		
		// did we skip something?
		if (!writer.writeWarnings().isEmpty()) {
			String details = "";
			for (String warning: writer.writeWarnings()) {
				details += warning + "\n";
			}
			application.showDialog("Workflow not fully saved", "Some parts of workflow structure are not supported by current workflow system and they were skipped. The rest of the workflow was successfully saved.", details, Severity.INFO, false);
		}
	}

	private void saveScript(File scriptFile, StringBuffer currentWorkflowScript) throws IOException {
		PrintWriter scriptOut = null;
		try {
			scriptOut = new PrintWriter(new OutputStreamWriter(new FileOutputStream(scriptFile)));
			scriptOut.print(currentWorkflowScript.toString());
			
		} finally {
			if (scriptOut != null) {
				scriptOut.close();
			}
		}
	}

	public File getScriptDirectory() {
		return scriptDirectory;
	}
}
