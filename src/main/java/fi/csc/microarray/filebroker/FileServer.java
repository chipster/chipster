package fi.csc.microarray.filebroker;

import java.io.File;
import java.io.IOException;
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
import fi.csc.microarray.messaging.message.NamiMessage;
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


    public FileServer() throws Exception {
    	
		// check file repository
		File fileRepository = new File(MicroarrayConfiguration.getValue("frontend", "fileServerPath"));
		if (!fileRepository.exists()) {
			boolean ok = fileRepository.mkdir();
			if (!ok) {
				throw new IOException("could not create file repository at " + fileRepository);
			}
		}

		// boot up file server
		EmbeddedJettyServer fileServer = new EmbeddedJettyServer(FILESERVER_CONTEXT_PATH);
		int port = FileBrokerConfig.getPort();
		fileServer.start(fileRepository.getPath(), "/", port);
		logger.info("fileserver is up and running [" + ApplicationConstants.NAMI_VERSION + "]");
		logger.info("[mem: " + MemUtil.getMemInfo() + "]");
		
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

    }


	public String getName() {
		return "filebroker";
	}


	public void onNamiMessage(NamiMessage msg) {
		// TODO Auto-generated method stub
		
	}
}
