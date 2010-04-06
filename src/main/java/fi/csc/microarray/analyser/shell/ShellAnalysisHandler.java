package fi.csc.microarray.analyser.shell;

import java.io.File;
import java.io.FileInputStream;
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

/**
 * Handles generic tools that are executed in a command line.
 * 
 * @author naktinis
 *
 */
public class ShellAnalysisHandler implements AnalysisHandler {

    private String descriptionDirectory;
    private String toolDirectory;
    
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
        ShellAnalysisJob analysisJob = new ShellAnalysisJob(toolDirectory, descriptionDirectory);
        analysisJob.construct(jobMessage, description, resultCallback);
        return analysisJob;
    }

    public AnalysisDescription handle(String descriptionFilename,
                                      Map<String, String> params)
            throws AnalysisException {
        
        // Generate analysis description
        AnalysisDescription ad = null;
        try {
            // FIXME path in configuration
            File sadlFile = new File(descriptionDirectory, descriptionFilename);
            byte[] bytes = new byte[(int) sadlFile.length()];
            new FileInputStream(sadlFile).read(bytes);
            String sadlString = new String(bytes);
            
            // Initiate description and set some basic values
            ad = new AnalysisDescriptionGenerator().generate(new ChipsterSADLParser().parse(sadlString), this);
            ad.setSADL(sadlString);
            
            // Command to be executed is stored in configuration file
            ad.setCommand(params.get("executable"));
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
