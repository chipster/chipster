package fi.csc.microarray.analyser.emboss;

import java.io.File;
import java.util.HashMap;

import org.apache.log4j.Logger;

import fi.csc.microarray.analyser.AnalysisDescription;
import fi.csc.microarray.analyser.AnalysisDescriptionGenerator;
import fi.csc.microarray.analyser.AnalysisException;
import fi.csc.microarray.analyser.AnalysisHandler;
import fi.csc.microarray.analyser.AnalysisJob;
import fi.csc.microarray.analyser.ResultCallback;
import fi.csc.microarray.analyser.bsh.BeanShellHandler;
import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.messaging.message.JobMessage;

/**
 * Responsible for creating AnalysisJobs and AnalysisDescriptions
 * for EMBOSS applications. Some parameters from configuration
 * should be known beforehand.
 * 
 * @author naktinis
 *
 */
public class EmbossAnalysisHandler implements AnalysisHandler {
    
    private String acdDirectory;
    private String toolDirectory;
    
    private static final Logger logger = Logger.getLogger(BeanShellHandler.class);

    /**
     * Constructor that takes configuration parameters
     * as input.
     * 
     * @param parameters
     */
    public EmbossAnalysisHandler(HashMap<String, String> parameters) {
        String externalToolPath = parameters.get("externalToolPath");
        String toolPath = parameters.get("toolPath");
        String descriptionPath = parameters.get("descriptionPath");
        
        toolDirectory = new File(externalToolPath, toolPath).getAbsolutePath();
        acdDirectory = new File(externalToolPath, descriptionPath).getAbsolutePath();
    }
    
    public AnalysisJob createAnalysisJob(JobMessage jobMessage, AnalysisDescription description,
                                         ResultCallback resultCallback) {
        EmbossAnalysisJob analysisJob = new EmbossAnalysisJob(toolDirectory, acdDirectory);
        analysisJob.construct(jobMessage, description, resultCallback);
        return analysisJob;
    }
    
    public AnalysisDescription handle(String acdFileName) throws AnalysisException {
        
        // Read ACD description
        File acdFile = new File(acdDirectory, acdFileName);
        ACDDescription acdDescription = new ACDDescription(acdFile);
        logger.debug("creating description from " + acdFile.getAbsolutePath());
        
        // Create description for computation server
        SADLDescription sadlDescription = new ACDToSADL(acdDescription).convert();
        AnalysisDescription description =
            new AnalysisDescriptionGenerator().generate(sadlDescription, this);
        
        // Fill description with Emboss-specific values
        description.setCommand("EMBOSS");
        description.setVVSADL(sadlDescription.toString());
        description.setSourceCode("Source code not available for EMBOSS tools. " + 
                                  "See: http://emboss.sourceforge.net/.");
        description.setSourceResourceName(acdFileName);
        description.setSourceResourceFullPath(acdFile.getAbsolutePath());
        
        return description;
    }
    
    public boolean isUptodate(AnalysisDescription description) {
        return true;
    }

    public boolean isDisabled() {
        return false;
    }
}
