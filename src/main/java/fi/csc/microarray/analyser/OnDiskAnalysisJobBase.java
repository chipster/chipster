package fi.csc.microarray.analyser;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.List;

import org.apache.log4j.Logger;

import fi.csc.microarray.analyser.ToolDescription.OutputDescription;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.filebroker.FileBrokerClient.FileBrokerArea;
import fi.csc.microarray.filebroker.NotEnoughDiskSpaceException;
import fi.csc.microarray.messaging.JobState;
import fi.csc.microarray.messaging.message.JobMessage;
import fi.csc.microarray.security.CryptoKey;
import fi.csc.microarray.util.Exceptions;
import fi.csc.microarray.util.Files;
import fi.csc.microarray.util.ToolUtils;

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
	public void construct(JobMessage inputMessage, ToolDescription analysis, ResultCallback resultHandler) {
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
		try {
			
			LinkedHashMap<String, String> nameMap = new LinkedHashMap<>();
			
			for (String fileName : inputMessage.getKeys()) {
				cancelCheck();

				// get url and output file
				String dataId = inputMessage.getId(fileName);
				File localFile = new File(jobWorkDir, fileName);
				
				// make local file available, by downloading, copying or symlinking
				resultHandler.getFileBrokerClient().getFile(dataId, new File(jobWorkDir, fileName));
				logger.debug("made available local file: " + localFile.getName() + " " + localFile.length());
				
				nameMap.put(fileName, inputMessage.getName(fileName));
			}
			
			ToolUtils.writeInputDescription(new File(jobWorkDir, "chipster-inputs.tsv"), nameMap);
			
		} catch (Exception e) {
			outputMessage.setErrorMessage("Transferring input data to computing service failed.");
			outputMessage.setOutputText(Exceptions.getStackTrace(e));
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
			
			// parse a file containing 
			String outputsFilename = "chipster-outputs.tsv";
			LinkedHashMap<String, String> nameMap = new LinkedHashMap<>();
			try {
				nameMap = ToolUtils.parseOutputDescription(new File(jobWorkDir, outputsFilename));
			} catch (IOException | MicroarrayException e) {
				logger.warn("couldn't parse " + outputsFilename);
				outputMessage.setErrorMessage("couldn't parse " + outputsFilename);
				outputMessage.setOutputText(Exceptions.getStackTrace(e));
                updateState(JobState.ERROR, "");
			}
			
			// add all described files to the result message
			for (File outputFile : describedFiles) {
	            // copy file to file broker
	            String dataId = CryptoKey.generateRandom();
	            try {
	                resultHandler.getFileBrokerClient().addFile(dataId, FileBrokerArea.CACHE, outputFile, null);
	                String nameInClient = nameMap.get(outputFile.getName());
	                // put dataId to result message
	                outputMessage.addPayload(outputFile.getName(), dataId, nameInClient);
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
	                
	            } catch (NotEnoughDiskSpaceException nedse) {
	            	logger.warn("not enough disk space for result file in filebroker");
	            	outputMessage.setErrorMessage("There was not enough disk space for the result file in the Chipster server. Please try again later.");
	            	updateState(JobState.FAILED_USER_ERROR, "not enough disk space for results");
	            }
	            
	            catch (Exception e) {
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
