package fi.csc.chipster.toolbox;

import java.io.IOException;
import java.util.Arrays;

import javax.jms.JMSException;
import javax.xml.parsers.ParserConfigurationException;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.xml.sax.SAXException;

import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.constants.ApplicationConstants;
import fi.csc.microarray.messaging.JMSMessagingEndpoint;
import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.MessagingListener;
import fi.csc.microarray.messaging.MessagingTopic;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.messaging.MonitoredNodeBase;
import fi.csc.microarray.messaging.Topics;
import fi.csc.microarray.messaging.message.ChipsterMessage;
import fi.csc.microarray.messaging.message.CommandMessage;
import fi.csc.microarray.messaging.message.ModuleDescriptionMessage;
import fi.csc.microarray.messaging.message.SourceMessage;
import fi.csc.microarray.service.KeepAliveShutdownHandler;
import fi.csc.microarray.service.ShutdownCallback;
import fi.csc.microarray.util.SystemMonitorUtil;

public class ToolboxServer extends MonitoredNodeBase implements MessagingListener, ShutdownCallback {

	private Logger logger = LogManager.getLogger();

	private MessagingEndpoint endpoint;
	private ToolboxService toolboxRestService;

	public ToolboxServer(String configURL) throws Exception {

		// get configs
		DirectoryLayout.initialiseServerLayout(Arrays.asList(new String[] { "toolbox" }), configURL);
		Configuration configuration = DirectoryLayout.getInstance().getConfiguration();

		// init logger
		logger.info("starting toolbox service...");

		String url = configuration.getString("messaging", "toolbox-url");
		String toolsBinPath = configuration.getString("toolbox", "tools-bin-path");

		this.toolboxRestService = new ToolboxService(url, toolsBinPath);
		this.toolboxRestService.startServerOldChipster();

		// initialize communications
		this.endpoint = new JMSMessagingEndpoint(this);
		MessagingTopic compTopic = endpoint.createTopic(Topics.Name.TOOLBOX_TOPIC, AccessMode.READ);
		compTopic.setListener(this);

		// create keep-alive thread and register shutdown hook
		KeepAliveShutdownHandler.init(this);

		logger.info("toolbox is up and running [" + ApplicationConstants.VERSION + "]");
		logger.info("[mem: " + SystemMonitorUtil.getMemInfo() + "]");
	}

	@Override
	public void onChipsterMessage(ChipsterMessage chipsterMessage) {

		// expect only command messages
		if (chipsterMessage instanceof CommandMessage) {
			CommandMessage commandMessage = (CommandMessage) chipsterMessage;

			// Request to send descriptions
			if (CommandMessage.COMMAND_DESCRIBE.equals(commandMessage.getCommand())) {
				try {
					handleGetDescriptions(commandMessage);
				} catch (Exception e) {
					logger.error("sending descriptions messages failed", e);
				}
				return;
			}

			// source code request
			else if (CommandMessage.COMMAND_GET_SOURCE.equals(commandMessage.getCommand())) {
				handleGetSourceCode(commandMessage);
				return;
			}

			// unknown command message
			else {
				logger.warn("unexpected command message: " + commandMessage.getCommand());
			}
		}

		// unknown message (something else than command message)
		else {
			logger.warn("expecting a command message, got something else: " + chipsterMessage.getMessageID() + " "
					+ chipsterMessage.getClass());
		}

	}

	@Override
	public String getName() {
		return "toolbox";
	}

	@Override
	public void shutdown() {
		logger.info("shutdown requested");

		// close toolbox rest
		try {
			this.toolboxRestService.close();
		} catch (Exception e) {
			logger.warn("closing toolbox rest service failed", e);
		}

		// close jms messaging endpoint
		try {
			this.endpoint.close();
		} catch (Exception e) {
			logger.warn("closing messaging endpoint failed", e);
		}

		logger.info("shutting down");
	}

	private void handleGetSourceCode(CommandMessage requestMessage) {
		String toolID = new String(requestMessage.getParameters().get(0));
		if (toolID == null || toolID.isEmpty()) {
			logger.warn("invalid source code request, tool id is: " + toolID);
			return;
		}

		ToolboxTool tool = toolboxRestService.getToolbox().getTool(toolID);

		if (tool != null) {
			// getSource() would return the original file without SADL replacements
			String source = tool.getSadlString() + tool.getCode();

			logger.info("sending source code for " + toolID);
			SourceMessage sourceMessage = new SourceMessage(source);
			if (sourceMessage != null) {
				sendReplyMessage(requestMessage, sourceMessage);
				return;
			}
		} else {
			logger.warn("tool " + toolID + " not found");
		}
	}

	private void handleGetDescriptions(CommandMessage requestMessage)
			throws IOException, SAXException, ParserConfigurationException {
		logger.info("sending all descriptions");

		for (ModuleDescriptionMessage descriptionMessage : toolboxRestService.getToolbox().getModuleDescriptions()) {
			descriptionMessage.setReplyTo(requestMessage.getReplyTo());
			logger.info("sending descriptions for module " + descriptionMessage.getModuleName());
			sendReplyMessage(requestMessage, descriptionMessage);

		}
	}

	/**
	 * Sends the message in new thread.
	 */
	private void sendReplyMessage(final ChipsterMessage original, final ChipsterMessage reply) {

		reply.setReplyTo(original.getReplyTo());

		new Thread(new Runnable() {
			public void run() {
				try {
					endpoint.replyToMessage(original, reply);
				} catch (JMSException e) {
					// Failing is ok, if some other comp has replied quicker and
					// the TempTopic has already been deleted
					// logger.error("Could not send message.", e);
				}
			}
		}).start();
	}

	public static void main(String[] args) throws Exception {
		new ToolboxServer(null);
	}
}
