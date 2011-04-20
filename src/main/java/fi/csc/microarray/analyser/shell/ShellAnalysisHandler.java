package fi.csc.microarray.analyser.shell;

import java.io.File;
import java.io.InputStream;
import java.util.HashMap;
import java.util.Map;

import org.apache.log4j.Logger;

import fi.csc.microarray.analyser.AnalysisDescription;
import fi.csc.microarray.analyser.AnalysisDescriptionGenerator;
import fi.csc.microarray.analyser.AnalysisException;
import fi.csc.microarray.analyser.AnalysisHandler;
import fi.csc.microarray.analyser.AnalysisJob;
import fi.csc.microarray.analyser.RepositoryModule;
import fi.csc.microarray.analyser.ResultCallback;
import fi.csc.microarray.messaging.message.JobMessage;
import fi.csc.microarray.module.chipster.ChipsterSADLParser;
import fi.csc.microarray.util.Files;

/**
 * Handles generic tools that are executed in a command line.
 * 
 * @author naktinis
 *
 */
public class ShellAnalysisHandler implements AnalysisHandler {
    
    private String DIRECTORY_NAME = "shell";

    private String toolPath;
    
    private static final Logger logger = Logger.getLogger(ShellAnalysisHandler.class);
    
    /**
     * Constructor that takes configuration parameters
     * as input.
     * 
     * @param parameters
     */
    public ShellAnalysisHandler(HashMap<String, String> parameters) {
        toolPath = parameters.get("toolPath");
    }
    
    public AnalysisJob createAnalysisJob(JobMessage jobMessage,
            AnalysisDescription description, ResultCallback resultCallback)
            throws AnalysisException {
        ShellAnalysisJob analysisJob = new ShellAnalysisJob();
        analysisJob.construct(jobMessage, description, resultCallback);
        return analysisJob;
    }

    public AnalysisDescription handle(RepositoryModule module, String descriptionFilename, Map<String, String> params)
            throws AnalysisException {
        
        // Generate analysis description
        AnalysisDescription ad = null;
        try {
    		File sadlFile = new File(module.getModuleDir(), toolPath + File.separator + descriptionFilename);

            String sadlString;
            if (sadlFile.exists()) {
                // Try opening a file using file system
                sadlString = Files.fileToString(sadlFile);
            } else {
                // Open file as a resource
                InputStream scriptSource = 
                    this.getClass().getResourceAsStream("/" + DIRECTORY_NAME + "/"
                                                        + descriptionFilename);
                sadlString = Files.inputStreamToString(scriptSource);
            }
            
            // Initiate description and set some basic values
            ad = new AnalysisDescriptionGenerator().generate(
            		new ChipsterSADLParser().parse(sadlString), this, module);
            ad.setSADL(sadlString);
            ad.setSourceCode(sadlString);
            
            // Command to be executed is stored in configuration file
            ad.setCommand(params.get("executable"));
            ad.setConfigParameters(params);
            
            // Log success
            logger.info("successfully loaded shell analysis description " + descriptionFilename);
        } catch (Exception e) {
            throw new AnalysisException(e);
        }
        return ad;
    }

    public boolean isDisabled() {
        return false;
    }

	/**
	 * Check if the source file has been modified since the 
	 * AnalysisDescription was created.
	 */
	public boolean isUptodate(AnalysisDescription description) {
		File scriptFile = description.getToolFile();
		return scriptFile.lastModified() <= description.getCreationTime();
	}
}
