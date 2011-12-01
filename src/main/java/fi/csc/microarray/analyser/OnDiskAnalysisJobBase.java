package fi.csc.microarray.analyser;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.net.URL;
import java.util.List;

import org.apache.log4j.Logger;

import fi.csc.microarray.analyser.AnalysisDescription.OutputDescription;
import fi.csc.microarray.messaging.JobState;
import fi.csc.microarray.messaging.message.JobMessage;
import fi.csc.microarray.util.Exceptions;
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
			outputMessage.setOutputText("Mitä perkelettä" + Exceptions.getStackTrace(e));
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
	    // update job state on the client side
		updateStateDetailToClient("transferring output data");
		cancelCheck();

		// pass output files to result message
		List<OutputDescription> outputFiles = analysis.getOutputFiles();
		for (OutputDescription fileDescription : outputFiles) {
			cancelCheck();
			
			// single file description can also describe several files
			File[] describedFiles;

			if (fileDescription.getFileName().isSpliced()) {
                // it is a set of files
			    String prefix = fileDescription.getFileName().getPrefix();
			    String postfix = fileDescription.getFileName().getPostfix();
			    String regex = prefix + ".*" + postfix;
			    describedFiles = Files.findFiles(jobWorkDir, regex);
			    
			    // if output is required there should be at least one
			    if (!fileDescription.isOptional() && describedFiles.length == 0) {
                    logger.error("required output file set not found");
                    outputMessage.setErrorMessage("Required output file set " +
                            fileDescription.getFileName().getID() + " is missing.");
                    updateState(JobState.ERROR, "");
                    return;
			    }
			} else {
			    // it is a single file
	            String outputName = fileDescription.getFileName().getID();
			    describedFiles = new File[] {new File(jobWorkDir, outputName)};
			}
			
			// add all described files to the result message
			for (File outputFile : describedFiles) {
	            // copy file to file broker
	            URL url;
	            try {
	                url = resultHandler.getFileBrokerClient().addFile(new FileInputStream(outputFile), null);
	                // put url to result message
	                outputMessage.addPayload(outputFile.getName(), url);
	                logger.debug("transferred output file: " + fileDescription.getFileName());

	            } catch (FileNotFoundException e) {
	                // required output file not found
	                if (!fileDescription.isOptional()) {
	                    logger.error("required output file not found", e);
	                    outputMessage.setErrorMessage("Required output file is missing.");
	                    outputMessage.setOutputText(Exceptions.getStackTrace(e));
	                    updateState(JobState.ERROR, "");
	                    return;
	                }
	                
	            } catch (Exception e) {
	                // TODO continue or return? also note the super.postExecute()
	                logger.error("could not put file to file broker", e);
	                outputMessage.setErrorMessage("Could not send output file.");
	                outputMessage.setOutputText(Exceptions.getStackTrace(e));
	                updateState(JobState.ERROR, "");
	                return;
	            }
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
