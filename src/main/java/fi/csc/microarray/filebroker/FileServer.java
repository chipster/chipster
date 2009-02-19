package fi.csc.microarray.filebroker;

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.Timer;

import org.apache.log4j.Logger;

import fi.csc.microarray.ApplicationConstants;
import fi.csc.microarray.MicroarrayConfiguration;
import fi.csc.microarray.manager.ManagerClient;
import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.MessagingListener;
import fi.csc.microarray.messaging.MessagingTopic;
import fi.csc.microarray.messaging.NodeBase;
import fi.csc.microarray.messaging.Topics;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.messaging.message.CommandMessage;
import fi.csc.microarray.messaging.message.NamiMessage;
import fi.csc.microarray.messaging.message.UrlMessage;
import fi.csc.microarray.util.FileCleanUpTimerTask;
import fi.csc.microarray.util.MemUtil;

public class FileServer extends NodeBase implements MessagingListener {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(FileServer.class);

	private MessagingEndpoint endpoint;
	private ManagerClient managerClient;
	private AuthorisedUrlRepository urlRepository;


    public FileServer() {

    	try {
    		// check file repository
    		File fileRepository = new File(MicroarrayConfiguration.getValue("frontend", "fileServerPath"));
    		if (!fileRepository.exists()) {
    			boolean ok = fileRepository.mkdir();
    			if (!ok) {
    				throw new IOException("could not create file repository at " + fileRepository);
    			}
    		}

    		// initialise url repository
    		String host = MicroarrayConfiguration.getValue("filebroker", "url");
    		int port = FileBrokerConfig.getPort();
    		this.urlRepository = new AuthorisedUrlRepository(host, port);

    		// boot up file server
    		JettyFileServer fileServer = new JettyFileServer(urlRepository);
    		fileServer.start(fileRepository.getPath(), port);

    		// start scheduler
    		int cutoff = 1000 * Integer.parseInt(MicroarrayConfiguration.getValue("frontend", "fileLifeTime"));
    		int cleanUpFrequency = 1000 * Integer.parseInt(MicroarrayConfiguration.getValue("frontend", "cleanUpFrequency"));
    		int checkFrequency = 1000 * 5;
    		Timer t = new Timer("frontend-scheduled-tasks", true);
    		t.schedule(new FileCleanUpTimerTask(fileRepository, cutoff), 0, cleanUpFrequency);
    		t.schedule(new JettyCheckTimerTask(fileServer), 0, checkFrequency);

    		// initialise messaging
    		this.endpoint = new MessagingEndpoint(this);
    		MessagingTopic urlRequestTopic = endpoint.createTopic(Topics.Name.AUTHORISED_URL_TOPIC, AccessMode.READ);
    		urlRequestTopic.setListener(this);
    		this.managerClient = new ManagerClient(endpoint); 

    		// all done
    		logger.info("fileserver is up and running [" + ApplicationConstants.NAMI_VERSION + "]");
    		logger.info("[mem: " + MemUtil.getMemInfo() + "]");

    	} catch (Exception e) {
    		logger.error(e, e);
    	}
    }


	public String getName() {
		return "filebroker";
	}


	public void onNamiMessage(NamiMessage msg) {
		try {

			if (msg instanceof CommandMessage && CommandMessage.COMMAND_URL_REQUEST.equals(((CommandMessage)msg).getCommand())) {
				URL url = urlRepository.createAuthorisedUrl();
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
}
