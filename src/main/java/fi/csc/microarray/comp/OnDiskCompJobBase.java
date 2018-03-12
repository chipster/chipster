package fi.csc.microarray.comp;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.UUID;

import org.apache.log4j.Logger;

import fi.csc.microarray.comp.ToolDescription.OutputDescription;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.filebroker.ChecksumException;
import fi.csc.microarray.filebroker.FileBrokerClient.FileBrokerArea;
import fi.csc.microarray.filebroker.FileBrokerException;
import fi.csc.microarray.filebroker.NotEnoughDiskSpaceException;
import fi.csc.microarray.messaging.JobState;
import fi.csc.microarray.messaging.message.GenericJobMessage;
import fi.csc.microarray.util.Exceptions;
import fi.csc.microarray.util.Files;
import fi.csc.microarray.util.ToolUtils;

/**
 * Provides functionality for transferring input files from file broker
 * to job work directory and output files from job work directory to file 
 * broker.
 *
 */
public abstract class OnDiskCompJobBase extends CompJob {

	private static final Logger logger = Logger.getLogger(OnDiskCompJobBase.class);

	private static final String JOB_DATA_DIR_NAME = "data";
	private static final String JOB_TOOLBOX_DIR_NAME = "toolbox";
	
	protected File jobDir;
	protected File jobDataDir;
	protected File jobToolboxDir;
	
	@Override
	public void construct(GenericJobMessage inputMessage, ToolDescription toolDescription, ResultCallback resultHandler) {
		super.construct(inputMessage, toolDescription, resultHandler);
		this.jobDir = new File(resultHandler.getWorkDir(), getId());
		this.jobDataDir = new File(this.jobDir, JOB_DATA_DIR_NAME);
		this.jobToolboxDir = new File(this.jobDir, JOB_TOOLBOX_DIR_NAME);
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

		updateState(JobState.RUNNING, "transferring input data");

		// create directories for the job
		if (!this.jobDir.mkdir()) {
			this.setErrorMessage("Creating job working directory failed.");
			updateState(JobState.ERROR);
			return;
		}

		// get input files and toolbox
		try {
			// input files
			getInputFiles();
			
			// toolbox
			if (!this.jobToolboxDir.mkdir()) {
				throw new IOException("Creating job toolbox dir failed.");
			}
			resultHandler.getToolboxClient().getToolboxModules(this.jobToolboxDir);

		} catch (Exception e) {
			this.setErrorMessage("Transferring input data and tools to computing service failed.");
			this.setOutputText(Exceptions.getStackTrace(e));
			updateState(JobState.ERROR);
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
		updateState(JobState.RUNNING, "transferring output data");
		cancelCheck();

		// pass output files to result message
		List<OutputDescription> outputFiles = toolDescription.getOutputFiles();
		for (OutputDescription fileDescription : outputFiles) {
			cancelCheck();
			
			// single file description can also describe several files
			File[] describedFiles;

			if (fileDescription.getFileName().isSpliced()) {
                // it is a set of files
			    String prefix = fileDescription.getFileName().getPrefix();
			    String postfix = fileDescription.getFileName().getPostfix();
			    String regex = prefix + ".*" + postfix;
			    describedFiles = Files.findFiles(jobDataDir, regex);
			    
			    // if output is required there should be at least one
			    if (!fileDescription.isOptional() && describedFiles.length == 0) {
                    logger.error("required output file set not found");
                    this.setErrorMessage("Required output file set " +
                            fileDescription.getFileName().getID() + " is missing.");
                    updateState(JobState.ERROR);
                    return;
			    }
			} else {
			    // it is a single file
	            String outputName = fileDescription.getFileName().getID();
			    describedFiles = new File[] {new File(jobDataDir, outputName)};
			}
			
			// parse a file containing file names for the client
			String outputsFilename = "chipster-outputs.tsv";
			LinkedHashMap<String, String> nameMap = new LinkedHashMap<>();
			try {
				nameMap = ToolUtils.parseOutputDescription(new File(jobDataDir, outputsFilename));
			} catch (IOException | MicroarrayException e) {
				logger.warn("couldn't parse " + outputsFilename);
				this.setErrorMessage("could not parse " + outputsFilename);
				this.setOutputText(Exceptions.getStackTrace(e));
                updateState(JobState.ERROR);
			}
			
			// add all described files to the result message
			for (File outputFile : describedFiles) {
	            // copy file to file broker
	            try {	            	
	            	String nameInClient = nameMap.get(outputFile.getName());
	            	String nameInSessionDb = nameInClient != null? nameInClient : outputFile.getName();
	                String dataId = resultHandler.getFileBrokerClient().addFile(UUID.fromString(inputMessage.getJobId()), inputMessage.getSessionId(), FileBrokerArea.CACHE, outputFile, null, nameInSessionDb);
	                // put dataId to result message
	                this.addOutputDataset(outputFile.getName(), dataId, nameInClient);
	                logger.debug("transferred output file: " + fileDescription.getFileName());

	            } catch (FileNotFoundException e) {
	                // required output file not found
	                if (!fileDescription.isOptional()) {
	                    logger.error("required output file not found", e);
	                    this.setErrorMessage("Required output file is missing.");
	                    this.appendOutputText(Exceptions.getStackTrace(e));
	                    updateState(JobState.ERROR);
	                    return;
	                }
	                
	            } catch (NotEnoughDiskSpaceException nedse) {
	            	logger.warn("not enough disk space for result file in filebroker");
	            	this.setErrorMessage("There was not enough disk space for the result file in the Chipster server. Please try again later.");
	            	updateState(JobState.FAILED_USER_ERROR, "not enough disk space for results");
	            }
	            
	            catch (Exception e) {
	                // TODO continue or return? also note the super.postExecute()
	                logger.error("could not put file to file broker", e);
	                this.setErrorMessage("Could not send output file.");
	                this.setOutputText(Exceptions.getStackTrace(e));
	                updateState(JobState.ERROR);
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
				Files.delTree(jobDir);
			}
		} catch (Exception e) {
			logger.error("Error when cleaning up job work dir.", e);
		} finally {
			super.cleanUp();
		}
	}


	private void getInputFiles()
			throws Exception, JobCancelledException, IOException, FileBrokerException, ChecksumException {
		LinkedHashMap<String, String> nameMap = new LinkedHashMap<>();
		
		if (!this.jobDataDir.mkdir()) {
			throw new IOException("Creating job data dir failed.");
		}
		
		for (String fileName : inputMessage.getKeys()) {
			cancelCheck();
	
			// get url and output file
			String dataId = inputMessage.getId(fileName);
			File localFile = new File(jobDataDir, fileName);
			
			// make local file available, by downloading, copying or symlinking
			resultHandler.getFileBrokerClient().getFile(inputMessage.getSessionId(), dataId, new File(jobDataDir, fileName));
			logger.debug("made available local file: " + localFile.getName() + " " + localFile.length());
			
			nameMap.put(fileName, inputMessage.getName(fileName));
		}
		
		ToolUtils.writeInputDescription(new File(jobDataDir, "chipster-inputs.tsv"), nameMap);
	
		inputMessage.preExecute(jobDataDir);
	}
	
}
