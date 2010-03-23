package fi.csc.microarray.analyser;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.net.URL;
import java.util.List;

import org.apache.log4j.Logger;

import fi.csc.microarray.messaging.JobState;
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
	 * @throws JobCancelledException 
	 * 
	 */
	@Override
	protected void preExecute() throws JobCancelledException {
		cancelCheck();
		super.preExecute();

		updateStateDetailToClient("transferring input data");

		// create working dir for the job
		if (!this.jobWorkDir.mkdir()) {
			outputMessage.setErrorMessage("Creating working directory failed.");
			updateState(JobState.ERROR, "");
			return;
		}

		// extract input files to work dir
		// TODO security check input file names
		try {
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
		} catch (Exception e) {
			outputMessage.setErrorMessage("Transferring input data to computing service failed.");
			outputMessage.setOutputText(e.toString());
			updateState(JobState.ERROR, "");
			return;
		}			
	}

	
	/**
	 * Copy output files from job work dir to file broker.
	 * 
	 */
	@Override
	protected void postExecute() throws JobCancelledException {
		updateStateDetailToClient("transferring output data");
		cancelCheck();

		List<String> outputFileNames = analysis.getOutputFiles();
		for (String fileName : outputFileNames) {
			cancelCheck();

		List<String> outputFileNames = analysis.getOutputFiles();
		for (String fileName : outputFileNames) {
			cancelCheck();

			// copy file to file broker
			File outputFile = new File(jobWorkDir, fileName);
			URL url;
				try {
					url = resultHandler.getFileBrokerClient().addFile(new FileInputStream(outputFile), null);
					// put url to result message
					outputMessage.addPayload(fileName, url);
					logger.debug("transferred output file: " + fileName);
						

			} catch (FileNotFoundException e) {
			    // FIXME need to deal missing required files
			    // Output file not found, it might have been optional.
			    // In future we might consider displaying an error message
			    // when a required output was not found.

			} catch (Exception e) {
					// TODO continue or return? also note the super.postExecute()
					logger.error("could not put file to file broker", e);
					outputMessage.setErrorMessage("Could not send output file.");
					outputMessage.setOutputText(e.toString());
					updateState(JobState.ERROR, "");
					return;
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
