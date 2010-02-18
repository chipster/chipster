package fi.csc.microarray.analyser.emboss;

import java.io.File;

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

public class EmbossAnalysisHandler implements AnalysisHandler {
    // FIXME temporary storage, use configuration file later
    private static final String scriptDirectory = "src/test/resources/";
    private static final Logger logger = Logger.getLogger(BeanShellHandler.class);

    
    public AnalysisJob createAnalysisJob(JobMessage jobMessage, AnalysisDescription description,
                                         ResultCallback resultCallback) {
        EmbossAnalysisJob analysisJob = new EmbossAnalysisJob();
        analysisJob.construct(jobMessage, description, resultCallback);
        return analysisJob;
    }
    
    public AnalysisDescription handle(String acdFileName) throws AnalysisException {
        
        // Read ACD description
        String acdPath = scriptDirectory + File.separator + acdFileName;
        File acdFile = new File(acdPath);
        ACDDescription acdDescription = new ACDDescription(acdFile);
        logger.debug("creating description from " + acdPath);
        
        // Create description for computation server
        SADLDescription sadlDescription = new ACDToSADL(acdDescription).convert();
        AnalysisDescription description =
            new AnalysisDescriptionGenerator().generate(sadlDescription, null);
        
        // Fill description with Emboss-specific values
        description.setCommand("EMBOSS");
        description.setSourceCode("Source code not available for EMBOSS tools. " + 
                                  "See: http://emboss.sourceforge.net/ ");
        description.setSourceResourceName(acdFileName);
        description.setSourceResourceFullPath(acdPath);
        
        return description;
    }
    
    /**
     * Check if the source file has been modified since the 
     * AnalysisDescription was created.
     */
    public boolean isUptodate(AnalysisDescription description) {
        return true;
    }

    public boolean isDisabled() {
        return false;
    }
}
