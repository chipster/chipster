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
import fi.csc.microarray.analyser.ResultCallback;
import fi.csc.microarray.description.SADLParser.ParseException;
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

    private String descriptionDirectory;
    
    private static final Logger logger = Logger.getLogger(ShellAnalysisHandler.class);
    
    /**
     * Constructor that takes configuration parameters
     * as input.
     * 
     * @param parameters
     */
    public ShellAnalysisHandler(HashMap<String, String> parameters) {
        descriptionDirectory = parameters.get("descriptionPath");
    }
    
    public AnalysisJob createAnalysisJob(JobMessage jobMessage,
            AnalysisDescription description, ResultCallback resultCallback)
            throws AnalysisException {
        ShellAnalysisJob analysisJob = new ShellAnalysisJob(description);
        analysisJob.construct(jobMessage, description, resultCallback);
        return analysisJob;
    }

    public AnalysisDescription handle(String descriptionFilename,
                                      Map<String, String> params)
            throws AnalysisException {
        
        // Generate analysis description
        AnalysisDescription ad = null;
        try {
            File sadlFile = new File(descriptionDirectory, descriptionFilename);
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
                    new ChipsterSADLParser().parse(sadlString), this);
            ad.setSADL(sadlString);
            
            // Command to be executed is stored in configuration file
            ad.setCommand(params.get("executable"));
            ad.setConfigParameters(params);
            
            // Log success
            logger.info("successfully loaded shell analysis description " + descriptionFilename);
        } catch (ParseException e) {
            throw new AnalysisException(e);
        } catch (Exception e) {
            e.printStackTrace();
        }
        return ad;
    }

    public boolean isDisabled() {
        return false;
    }

    public boolean isUptodate(AnalysisDescription description) {
        return true;
    }
}
