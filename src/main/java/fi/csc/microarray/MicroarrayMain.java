package fi.csc.microarray;


import java.io.FileInputStream;

import fi.csc.microarray.analyser.AnalyserServer;
import fi.csc.microarray.analyser.r.VVSADLTool;
import fi.csc.microarray.auth.Authenticator;
import fi.csc.microarray.client.SwingClientApplication;
import fi.csc.microarray.filebroker.FileServer;
import fi.csc.microarray.manager.Manager;
import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.NodeBase;
import fi.csc.microarray.messaging.Topics;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.module.chipster.ChipsterVVSADLParser.Validator;
import fi.csc.microarray.util.CommandLineParser;
import fi.csc.microarray.util.CommandLineParser.CommandLineException;
import fi.csc.microarray.util.config.ConfigurationLoader.OldConfigurationFormatException;
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
			cmdParser.addParameter("analyser", false, false, null, "start analyser");
			cmdParser.addParameter("webstart", false, false, null, "start webstart service");
			cmdParser.addParameter("manager", false, false, null, "start manager service");
			cmdParser.addParameter("tests", false, false, null, "run tests");
			cmdParser.addParameter("nagios-check", false, false, null, "do nagios-compatitible system availability check");			
			cmdParser.addParameter("system-status", false, false, null, "query and print system status");
			cmdParser.addParameter("broker-check", false, false, null, "check broker availability");
			cmdParser.addParameter("rcheck", false, true, null, "check R script syntax");
			cmdParser.addParameter("-override", false, true, null, "comma separated list of configuration overrides (e.g. \"-override ModuleA/EntryA=val1,ModuleA/EntryB=val1;val2\")");
			cmdParser.addParameter("-homework", false, false, null, "store work files in home instead of workdir");
			cmdParser.addParameter("-required-analyser-count", false, true, "1", "store work files in home instead of workdir");
			cmdParser.addParameter("-repository-url", false, true, null, "url for the fileserver repository for rest test");
			cmdParser.addParameter("-testfile", false, true, null, "testfile for rest test");
			cmdParser.addParameter("-threads", false, true, "1", "thread count for rest test");
			cmdParser.addParameter("-repeat", false, true, "1", "repeat count for rest test");
			cmdParser.addParameter("-put", false, false, null, "use put for rest test");
			cmdParser.addParameter("-get", false, false, null, "use get count for rest test");
			
			
			// parse commandline
			cmdParser.parse(args);
			
			// give help, if needed
			if (cmdParser.userAskedHelp()) {
				System.out.println("Chipster " + ApplicationConstants.NAMI_VERSION);
				System.out.println("Parameters:");
				System.out.println(cmdParser.getDescription());
				System.exit(0);
			}

			// load configuration, when needed
			boolean needConfig = !(cmdParser.hasValue("rcheck")); 
			if (needConfig) {
				String overrides = cmdParser.getValue("-override");
				boolean useHomeAsWorkDir = cmdParser.hasValue("-homework");
				boolean defaultsUsed;
				try {
					defaultsUsed = MicroarrayConfiguration.loadConfiguration(overrides, useHomeAsWorkDir);
				} catch (OldConfigurationFormatException e) {
					if (!cmdParser.hasValue("authenticator") && !cmdParser.hasValue("analyser") && !cmdParser.hasValue("fileserver")) {
						SwingClientApplication.reportOldConfigurationFormatException(e);
						return;
					} else {
						throw e;
					}
					
				}
				if (defaultsUsed && !cmdParser.hasValue("nagios-check")) {
					System.out.println("No configuration found, defaults used and new configuration file created.");
				}
			}

			// start application
			if (cmdParser.hasValue("authenticator")) {
				new Authenticator();
				
			} else if (cmdParser.hasValue("analyser")) {
				new AnalyserServer();

			} else if (cmdParser.hasValue("fileserver")) {
				new FileServer();
			
			} else if (cmdParser.hasValue("webstart")) {
				new WebstartJettyServer().start();
			
			} else if (cmdParser.hasValue("manager")) {
				new Manager();

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
				
//			} else if (cmdParser.hasValue("tests")) {
//				junit.textui.TestRunner.run(JmsTests.suite());
				
			} else if (cmdParser.hasValue("rcheck")) {
				boolean fails = false;
				try {					
					VVSADLTool.ParsedRScript res = new VVSADLTool().parseRScript(new FileInputStream(cmdParser.getValue("rcheck")));
					new Validator().validate(cmdParser.getValue("rcheck"), res.VVSADL);
				} catch (Exception e) {
					System.out.println(e.getMessage());
					fails = true;
				}
				System.out.println("parse succeeded: " + !fails);
				
			} else {
				SwingClientApplication.start();				
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
