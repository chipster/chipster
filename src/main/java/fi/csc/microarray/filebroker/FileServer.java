package fi.csc.microarray.filebroker;

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.Timer;

import org.apache.log4j.Logger;

import fi.csc.microarray.ApplicationConstants;
import fi.csc.microarray.MicroarrayConfiguration;
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

    private static final String FILESERVER_CONTEXT_PATH = "/fileserver";

	private MessagingEndpoint endpoint;
	private MessagingTopic managerTopic;
	private AuthorisedUrlRepository urlRepository;


    public FileServer() throws Exception {
    	
		// check file repository
		File fileRepository = new File(MicroarrayConfiguration.getValue("frontend", "fileServerPath"));
		if (!fileRepository.exists()) {
			boolean ok = fileRepository.mkdir();
			if (!ok) {
				throw new IOException("could not create file repository at " + fileRepository);
			}
		}

    	// initialise url repository
		URL rootUrl = new URL(MicroarrayConfiguration.getValue("filebroker", "url"));
    	this.urlRepository = new AuthorisedUrlRepository(rootUrl);

		// boot up file server
		JettyFileServer fileServer = new JettyFileServer(FILESERVER_CONTEXT_PATH, urlRepository, rootUrl);
		int port = FileBrokerConfig.getPort();
		fileServer.start(fileRepository.getPath(), "/", port);

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
		managerTopic = endpoint.createTopic(Topics.Name.MANAGER_TOPIC, AccessMode.WRITE);

		// all done
		logger.info("fileserver is up and running [" + ApplicationConstants.NAMI_VERSION + "]");
		logger.info("[mem: " + MemUtil.getMemInfo() + "]");
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
				
			} else {
				logger.error("message " + msg.getMessageID() + " not understood");
			}
			
		} catch (Exception e) {
			logger.error(e, e);
		}
	}
}
