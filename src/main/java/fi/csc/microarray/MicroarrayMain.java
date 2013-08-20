package fi.csc.microarray;


import java.io.FileInputStream;
import java.io.IOException;

import javax.swing.SwingUtilities;

import fi.csc.microarray.analyser.AnalyserServer;
import fi.csc.microarray.analyser.SADLTool;
import fi.csc.microarray.auth.Authenticator;
import fi.csc.microarray.client.SwingClientApplication;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.constants.ApplicationConstants;
import fi.csc.microarray.filebroker.FileServer;
import fi.csc.microarray.manager.Manager;
import fi.csc.microarray.messaging.AdminAPI;
import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.NodeBase;
import fi.csc.microarray.messaging.Topics;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
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
			cmdParser.addParameter("standalone", false, false, null, "start standalone client");
			cmdParser.addParameter("authenticator", false, false, null, "start authenticator");
			cmdParser.addParameter("fileserver", false, false, null, "start fileserver");
			cmdParser.addParameter("analyser", false, false, null, "start analyser");
			cmdParser.addParameter("webstart", false, false, null, "start webstart service");
			cmdParser.addParameter("manager", false, false, null, "start manager service");
			cmdParser.addParameter("tests", false, false, null, "run tests");
			cmdParser.addParameter("nagios-check", false, false, null, "do nagios-compatitible system availability check");			
			cmdParser.addParameter("system-status", false, false, null, "query and print system status");
			cmdParser.addParameter("broker-check", false, false, null, "check broker availability");
			cmdParser.addParameter("rcheck", false, true, null, "check R script syntax");
			cmdParser.addParameter("-config", false, true, null, "configuration file URL (chipster-config.xml)");
			cmdParser.addParameter("-required-analyser-count", false, true, "1", "required comp service count for nagios check");
            cmdParser.addParameter("-module", false, true, "fi.csc.microarray.module.chipster.MicroarrayModule", "client module (e.g. microarray-module)");
			
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
				
			} else if (cmdParser.hasValue("analyser")) {
				new AnalyserServer(configURL);

			} else if (cmdParser.hasValue("fileserver")) {
				new FileServer(configURL);
			
			} else if (cmdParser.hasValue("webstart")) {
				new WebstartJettyServer().start();
			
			} else if (cmdParser.hasValue("manager")) {
				new Manager(configURL);

			} else if (cmdParser.hasValue("nagios-check") || cmdParser.hasValue("system-status")) {
				
				// query status
				int requiredAnalyserCount = Integer.parseInt(cmdParser.getValue("-required-analyser-count"));
				boolean ok;
				String error = "";				
				String status = "";
				try {
					NodeBase nodeSupport = new NodeBase() {
						public String getName() {
							return "nagios-check";
						}
					};
					DirectoryLayout.initialiseSimpleLayout(configURL).getConfiguration();       			    
					MessagingEndpoint endpoint = new MessagingEndpoint(nodeSupport);
					AdminAPI api = new AdminAPI(endpoint.createTopic(Topics.Name.ADMIN_TOPIC, AccessMode.READ_WRITE), null);
					api.setRequiredCountFor("analyser", requiredAnalyserCount);
					boolean fastCheck = cmdParser.hasValue("nagios-check");
					ok = api.areAllServicesUp(fastCheck);
					error = api.getErrorStatus();
					status = api.generateStatusReport();
					endpoint.close();
					
				} catch (Exception e) {
					ok = false;
					error = e.getMessage();
				}
				
				// print results
				if (cmdParser.hasValue("nagios-check")) {
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
				
			}  else if (cmdParser.hasValue("broker-check")) {
				
				String error = "";
				
				try {
					NodeBase nodeSupport = new NodeBase() {
						public String getName() {
							return "nagios-check";
						}
					};
					MessagingEndpoint endpoint = new MessagingEndpoint(nodeSupport);
					endpoint.close();
					
				} catch (Exception e) {
					System.out.println("BROKER NOT AVAILABLE: " + error);
					System.exit(1);
				}
				System.out.println("broker available");
				System.exit(0);
				
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

			} else if (cmdParser.hasValue("standalone")) {
				
				final String module = cmdParser.getValue("-module");
				
				SwingUtilities.invokeLater(new Runnable() {
					@Override
					public void run() {						
						try {
							SwingClientApplication.startStandalone(module);
						} catch (IOException e) {
							e.printStackTrace();
							System.exit(0);
						}		
					}
				});
				
			} else {

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
