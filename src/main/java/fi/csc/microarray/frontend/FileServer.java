package fi.csc.microarray.frontend;

import java.io.File;
import java.io.IOException;
import java.util.Timer;

import org.apache.log4j.Logger;

import fi.csc.microarray.ApplicationConstants;
import fi.csc.microarray.MicroarrayConfiguration;
import fi.csc.microarray.util.FileCleanUpTimerTask;
import fi.csc.microarray.util.MemUtil;

public class FileServer {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(FileServer.class);

    private static final String FILESERVER_CONTEXT_PATH = "/fileserver";


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

    }
}
