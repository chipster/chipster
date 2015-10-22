package fi.csc.microarray;


import java.io.FileInputStream;
import java.io.IOException;

import javax.swing.SwingUtilities;

import fi.csc.chipster.toolbox.SADLTool;
import fi.csc.microarray.auth.Authenticator;
import fi.csc.microarray.client.SwingClientApplication;
import fi.csc.microarray.comp.CompServer;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.constants.ApplicationConstants;
import fi.csc.microarray.filebroker.FileServer;
import fi.csc.microarray.jobmanager.JobManager;
import fi.csc.microarray.manager.Manager;
import fi.csc.microarray.messaging.JMSMessagingEndpoint;
import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.messaging.NodeBase;
import fi.csc.microarray.messaging.Topics;
import fi.csc.microarray.messaging.admin.AdminAPI;
import fi.csc.microarray.module.chipster.ChipsterSADLParser.Validator;
import fi.csc.microarray.util.CommandLineParser;
import fi.csc.microarray.util.CommandLineParser.CommandLineException;
import fi.csc.microarray.webstart.WebstartJettyServer;

/**
 * The main program of Chipster system, actually just a loader for the 
 * actual components.
 * 
 * @author Aleksi Kallio, Taavi Hupponen
 */
public class MicroarrayMain {

	public static void main(String[] args) {

		try {

			// create args descriptions
			CommandLineParser cmdParser = new CommandLineParser();
			cmdParser.addParameter("client", false, false, null, "start client (default)");
			cmdParser.addParameter("authenticator", false, false, null, "start authenticator");
			cmdParser.addParameter("fileserver", false, false, null, "start fileserver");
			cmdParser.addParameter("comp", false, false, null, "start comp");
			cmdParser.addParameter("webstart", false, false, null, "start webstart service");
			cmdParser.addParameter("manager", false, false, null, "start manager service");
			cmdParser.addParameter("jobmanager", false, false, null, "start jobmanager service");
			cmdParser.addParameter("ping", false, false, null, "query and print system status");
			cmdParser.addParameter("ping-nagios", false, false, null, "query and print system status in nagios compatible format");			
			cmdParser.addParameter("rcheck", false, true, null, "check R script syntax");
			cmdParser.addParameter("-config", false, true, null, "configuration file URL (chipster-config.xml)");
            cmdParser.addParameter("-module", false, true, null, "client module (e.g. fi.csc.microarray.module.chipster.MicroarrayModule)");
			
			// parse commandline
			cmdParser.parse(args);
			
			// configuration file path
			String configURL = cmdParser.getValue("-config");
			
			// give help, if needed
			if (cmdParser.userAskedHelp()) {
				System.out.println("Chipster " + ApplicationConstants.VERSION);
				System.out.println("Parameters:");
				System.out.println(cmdParser.getDescription());
				System.exit(0);
			}

			// start application		
			
			if (cmdParser.hasValue("authenticator")) {
				new Authenticator(configURL);
				
			} else if (cmdParser.hasValue("comp")) {
				new CompServer(configURL);

			} else if (cmdParser.hasValue("fileserver")) {
				new FileServer(configURL, null);
			
			} else if (cmdParser.hasValue("webstart")) {
				new WebstartJettyServer().start();
			
			} else if (cmdParser.hasValue("manager")) {
				new Manager(configURL);

			} else if (cmdParser.hasValue("jobmanager")) {
				new JobManager(configURL);

			} else if (cmdParser.hasValue("ping") || cmdParser.hasValue("ping-nagios")) {
				
				// query status
				boolean ok;
				String error = "";				
				String status = "";
				try {
					NodeBase nodeSupport = new NodeBase() {
						public String getName() {
							return "ping";
						}
					};
					DirectoryLayout.initialiseSimpleLayout(configURL).getConfiguration();       			    
					MessagingEndpoint endpoint = new JMSMessagingEndpoint(nodeSupport);
					AdminAPI api = new AdminAPI(endpoint.createTopic(Topics.Name.ADMIN_TOPIC, AccessMode.READ_WRITE), null);
					boolean fastCheck = cmdParser.hasValue("ping-nagios");
					ok = api.areAllServicesUp(fastCheck);
					error = api.getErrorStatus();
					status = api.generateStatusReport();
					endpoint.close();
					
				} catch (Exception e) {
					ok = false;
					error = e.getMessage();
				}
				
				// print results
				if (cmdParser.hasValue("ping-nagios")) {
					if (ok) {
						System.out.println("CHIPSTER OK");
						System.exit(0);
					} else {
						System.out.println("CHIPSTER FAILED: " + error);
						System.exit(2);
					}
				} else {
					if (ok) {
						System.out.println("Chipster OK.");
					} else {
						System.out.println("Chipster failed: " + error);
					}
					System.out.println(status);
				}
				
			} else if (cmdParser.hasValue("rcheck")) {
				boolean fails = false;
				try {					
					SADLTool.ParsedScript res = new SADLTool("#").parseScript(new FileInputStream(cmdParser.getValue("rcheck")));
					new Validator().validate(cmdParser.getValue("rcheck"), res.SADL);
				} catch (Exception e) {
					System.out.println(e.getMessage());
					fails = true;
				}
				System.out.println("parse succeeded: " + !fails);
				
			} else {
				
				// assume client by default

				final String module = cmdParser.getValue("-module");
				final String config = configURL;								
				
				SwingUtilities.invokeLater(new Runnable() {
					@Override
					public void run() {						
						try {
							SwingClientApplication.start(config, module);		
						} catch (IOException e) {
							e.printStackTrace();
							System.exit(0);
						}		
					}
				});
			}
			
		} catch (CommandLineException e) {
			System.out.println("Illegal parameters");
			System.out.println("  " + e.getMessage());
			System.exit(1);
			
		} catch (Throwable t) {
			t.printStackTrace();
			System.exit(1);
		}		
	}
}
