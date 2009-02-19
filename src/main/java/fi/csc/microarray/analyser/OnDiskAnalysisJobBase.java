package fi.csc.microarray.analyser;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.net.URL;
import java.util.List;

import org.apache.log4j.Logger;

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

			logger.debug("Input file created: " + outputFile.getName() + " " + outputFile.length());
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
			// copy file to file broker
			File outputFile = new File(jobWorkDir, fileName);
			URL url = resultHandler.getFileBrokerClient().addFile(new FileInputStream(outputFile), null);
			outputMessage.addPayload(fileName, url);
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
