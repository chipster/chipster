package fi.csc.microarray.analyser;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.net.URL;
import java.util.List;

import org.apache.log4j.Logger;

import fi.csc.microarray.messaging.message.JobMessage;
import fi.csc.microarray.util.Files;
import fi.csc.microarray.util.IOUtils;

/**
 * Provides functionality for transferring input files from file broker
 * to job work directory and output files from job work directory to file 
 * broker.
 *
 */
public abstract class OnDiskAnalysisJobBase extends AnalysisJob {

	private static final Logger logger = Logger.getLogger(OnDiskAnalysisJobBase.class);

	protected File jobWorkDir;

	@Override
	public void construct(JobMessage inputMessage, AnalysisDescription analysis, ResultCallback resultHandler) {
		super.construct(inputMessage, analysis, resultHandler);
		this.jobWorkDir = new File(resultHandler.getWorkDir(), getId());
	}


	/**
	 * Copy input files from file broker to job work directory.
	 * 
	 */
	@Override
	protected void preExecute() throws Exception {
		cancelCheck();
		super.preExecute();

		updateStateDetail("transferring input data", true);

		// create working dir for the job
		if (!this.jobWorkDir.mkdir()) {
			throw new IOException("Could not create work dir: " + jobWorkDir.toString());
		}

		// extract input files to work dir
		// TODO security check input file names
		for (String fileName : inputMessage.payloadNames()) {
			cancelCheck();

			// get url
			URL inputUrl = inputMessage.getPayload(fileName);

			// get stream
			File outputFile;
			BufferedInputStream inputStream = null;
			BufferedOutputStream fileStream = null;
			try {
				inputStream = new BufferedInputStream(resultHandler.getFileBrokerClient().getFile(inputUrl));

				// copy to file
				outputFile = new File(jobWorkDir, fileName);
				fileStream = new BufferedOutputStream(new FileOutputStream(outputFile));
				IOUtils.copy(inputStream, fileStream);
			} finally {
				IOUtils.closeIfPossible(inputStream);
				IOUtils.closeIfPossible(fileStream);
			}

			logger.debug("created input file: " + outputFile.getName() + " " + outputFile.length());
		}
	}

	
	/**
	 * Copy output files from job work dir to file broker.
	 * 
	 */
	@Override
	protected void postExecute() throws Exception {
		updateStateDetail("transferring output data", true);
		cancelCheck();

		List<String> outputFileNames = analysis.getOutputFiles();
		for (String fileName : outputFileNames) {
			cancelCheck();

			try {
    			// Copy file to file broker
    			File outputFile = new File(jobWorkDir, fileName);
    			URL url = resultHandler.getFileBrokerClient().
    			              addFile(new FileInputStream(outputFile), null);
    
    			// Put url to result message
    			outputMessage.addPayload(fileName, url);
    			logger.debug("transferred output file: " + fileName);
			} catch (FileNotFoundException e) {
			    // Output file not found, it might have been optional.
			    // In future we might consider displaying an error message
			    // when a required output was not found.
			}
		}
		super.postExecute();
	}


	/**
	 * Clear job working directory.
	 * 
	 */
	@Override
	protected void cleanUp() {
		try {
			// sweep job working directory
			if (resultHandler.shouldSweepWorkDir()) {
				Files.delTree(jobWorkDir);
			}
		} catch (Exception e) {
			logger.error("Error when cleaning up job work dir.", e);
		} finally {
			super.cleanUp();
		}
	}
}
