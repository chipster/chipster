package fi.csc.microarray.analyser;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.List;

import org.apache.log4j.Logger;
import org.mortbay.util.IO;

import fi.csc.microarray.MicroarrayConfiguration;
import fi.csc.microarray.messaging.message.JobMessage;
import fi.csc.microarray.util.Files;
import fi.csc.microarray.util.IOUtils;

public abstract class OnDiskAnalysisJobBase extends AnalysisJob {

	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger
			.getLogger(OnDiskAnalysisJobBase.class);
	// FIXME this is not used
	public static final int TIMEOUT = Integer.parseInt(MicroarrayConfiguration
			.getValue("analyser", "timeout_sec")) * 1000;

	protected File jobWorkDir;

	@Override
	public void construct(JobMessage inputMessage,
			AnalysisDescription analysis, ResultCallback resultHandler) {
		super.construct(inputMessage, analysis, resultHandler);
		this.jobWorkDir = new File(resultHandler.getWorkDir(), getId());
	}

	@Override
	protected void preExecute() throws Exception {
		cancelCheck();
		super.preExecute();

		updateStateDetail("transferring input data", true);

		// create working dir
		if (!this.jobWorkDir.mkdir()) {
			throw new IOException("Could not create work dir: "
					+ jobWorkDir.toString());
		}

		// extract input files to work dir
		// TODO security check input file names
		logger.debug("Starting to write input files.");
		for (String fileName : inputMessage.payloadNames()) {
			cancelCheck();

			logger.debug("Getting input file: " + fileName);

			// check that payload exists on the file server
			URL inputURL = inputMessage.getPayloadURL(fileName);

			// make sure that the payload actually is available, as
			// jetty is sometimes a bit slow to write the payload to disk
			int maxWaitTime = 5120;
			int waitTime = 10;
			HttpURLConnection connection = null;
			try {
				connection  = (HttpURLConnection) inputURL.openConnection();
				int responseCode = connection.getResponseCode();
				while (responseCode != HttpURLConnection.HTTP_OK && waitTime <= maxWaitTime) {
					logger.info("Waiting for payload to become available.");
					waitTime = waitTime*2;
					try {
						Thread.sleep(waitTime);
					} catch (InterruptedException e) {
						logger.error("Interrupted while waiting for payload to become available.");
					}
					responseCode = connection.getResponseCode();
				}

				if (responseCode != HttpURLConnection.HTTP_OK) {
					throw new RetryException();
				}
			} finally {
				IOUtils.disconnectIfPossible(connection);
			}

				
			// payload exists on the server side fetch it
			File outputFile = new File(jobWorkDir, fileName);
			InputStream input = new BufferedInputStream(inputMessage
					.getPayload(fileName));

			OutputStream output = new BufferedOutputStream(
					new FileOutputStream(outputFile));
			IO.copy(input, output);
			input.close();
			output.flush();
			output.close();
			logger.debug("Input file created: " + outputFile.getName() + " "
					+ outputFile.length());
		}
		logger.debug("input files written");
	}

	@Override
	protected void postExecute() throws Exception {
		updateStateDetail("transferring output data", true);
		cancelCheck();

		List<String> outputFileNames = analysis.getOutputFiles();

		logger.debug("reading outputs of " + analysis.getFullName()
				+ ", total of " + outputFileNames.size() + " outputfiles");
		for (String fileName : outputFileNames) {
			cancelCheck();
			File outputFile = new File(jobWorkDir, fileName);
			outputMessage.addPayload(fileName, new BufferedInputStream(
					new FileInputStream(outputFile)));
			logger.debug("Added output file " + fileName);
		}

		super.postExecute();
	}

	protected void cleanUp() {
		try {
			// sweep working directory
			if (resultHandler.shouldSweepWorkDir()) {
				Files.delTree(jobWorkDir);
			}
		} catch (Exception e) {
			logger.error("Error in cleanUp.", e);
		} finally {
			super.cleanUp();
		}
	}

}
