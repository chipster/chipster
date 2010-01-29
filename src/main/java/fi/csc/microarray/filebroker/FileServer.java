package fi.csc.microarray.filebroker;

import java.io.File;
import java.net.URL;
import java.util.Arrays;
import java.util.Timer;

import javax.jms.JMSException;

import org.apache.log4j.Logger;

import fi.csc.microarray.ApplicationConstants;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.manager.ManagerClient;
import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.MessagingListener;
import fi.csc.microarray.messaging.MessagingTopic;
import fi.csc.microarray.messaging.NodeBase;
import fi.csc.microarray.messaging.Topics;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.messaging.message.CommandMessage;
import fi.csc.microarray.messaging.message.NamiMessage;
import fi.csc.microarray.messaging.message.ParameterMessage;
import fi.csc.microarray.messaging.message.UrlMessage;
import fi.csc.microarray.service.KeepAliveShutdownHandler;
import fi.csc.microarray.service.ShutdownCallback;
import fi.csc.microarray.util.FileCleanUpTimerTask;
import fi.csc.microarray.util.MemUtil;

public class FileServer extends NodeBase implements MessagingListener, ShutdownCallback {
	/**
	 * Logger for this class
	 */
	private static Logger logger;

	private MessagingEndpoint endpoint;	
	private ManagerClient managerClient;
	private AuthorisedUrlRepository urlRepository;


	public static void main(String[] args) {
		// we should be able to specify alternative user dir for testing... and replace maybe that previous hack
		//DirectoryLayout.getInstance(new File("chipster-userdir-fileserver")).getConfiguration();
		
		DirectoryLayout.getInstance().getConfiguration();
		new FileServer();
	}
	
    public FileServer() {

    	try {
    		// initialise dir and logging
    		DirectoryLayout.initialiseServerLayout(Arrays.asList(new String[] {"frontend", "filebroker"}));
    		Configuration configuration = DirectoryLayout.getInstance().getConfiguration();
    		logger = Logger.getLogger(FileServer.class);

    		// initialise url repository
    		File fileRepository = DirectoryLayout.getInstance().getFileRoot();
    		String host = configuration.getString("filebroker", "url");
    		int port = configuration.getInt("filebroker", "port");
    		
    		this.urlRepository = new AuthorisedUrlRepository(host, port);

    		// boot up file server
    		JettyFileServer fileServer = new JettyFileServer(urlRepository);
    		fileServer.start(fileRepository.getPath(), port);

    		// start scheduler
    		int cutoff = 1000 * configuration.getInt("filebroker", "file-life-time");
    		int cleanUpFrequency = 1000 * configuration.getInt("filebroker", "clean-up-frequency");
    		int checkFrequency = 1000 * 5;
    		Timer t = new Timer("frontend-scheduled-tasks", true);
    		t.schedule(new FileCleanUpTimerTask(fileRepository, cutoff), 0, cleanUpFrequency);
    		t.schedule(new JettyCheckTimerTask(fileServer), 0, checkFrequency);

    		// initialise messaging
    		this.endpoint = new MessagingEndpoint(this);
    		MessagingTopic urlRequestTopic = endpoint.createTopic(Topics.Name.AUTHORISED_URL_TOPIC, AccessMode.READ);
    		urlRequestTopic.setListener(this);
    		this.managerClient = new ManagerClient(endpoint); 

    		// create keep-alive thread and register shutdown hook
    		KeepAliveShutdownHandler.init(this);

    		
    		logger.info("fileserver is up and running [" + ApplicationConstants.NAMI_VERSION + "]");
    		logger.info("[mem: " + MemUtil.getMemInfo() + "]");

    	} catch (Exception e) {
    		e.printStackTrace();
    		logger.error(e, e);
    	}
    }


	public String getName() {
		return "filebroker";
	}


	public void onNamiMessage(NamiMessage msg) {
		try {

			if (msg instanceof CommandMessage && CommandMessage.COMMAND_URL_REQUEST.equals(((CommandMessage)msg).getCommand())) {
				CommandMessage requestMessage = (CommandMessage) msg;
				boolean useCompression = requestMessage.getParameters().contains(ParameterMessage.PARAMETER_USE_COMPRESSION);
				URL url = urlRepository.createAuthorisedUrl(useCompression);
				UrlMessage reply = new UrlMessage(url);
				endpoint.replyToMessage(msg, reply);
				managerClient.urlRequest(msg.getUsername(), url);
			} else {
				logger.error("message " + msg.getMessageID() + " not understood");
			}
			
		} catch (Exception e) {
			logger.error(e, e);
		}
	}

	public void shutdown() {
		logger.info("shutdown requested");

		// close messaging endpoint
		try {
			this.endpoint.close();
		} catch (JMSException e) {
			logger.error("closing messaging endpoint failed", e);
		}

		logger.info("shutting down");
	}


}
