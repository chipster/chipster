package fi.csc.microarray.analyser.emboss;

import java.io.File;
import java.util.HashMap;
import java.util.LinkedList;

import org.apache.log4j.Logger;

import fi.csc.microarray.analyser.OnDiskAnalysisJobBase;
import fi.csc.microarray.analyser.ResultCallback;
import fi.csc.microarray.analyser.AnalysisDescription.ParameterDescription;
import fi.csc.microarray.messaging.JobState;
import fi.csc.microarray.messaging.message.ResultMessage;

/**
 * Runs EMBOSS applications.
 * 
 * @author naktinis, akallio
 *
 */
public class EmbossAnalysisJob extends OnDiskAnalysisJobBase {
    
    private String toolDirectory;
    private String descriptionDirectory;
    
    // EMBOSS logger
    static final Logger logger = Logger.getLogger(EmbossAnalysisJob.class);
    
    LinkedList<EmbossQualifier> qualifiers;
    LinkedList<String> inputParameters;
    ACDDescription acdDescription;
    
    @Override
    protected void cancelRequested() {
        // TODO Auto-generated method stub

    }
    
    public EmbossAnalysisJob(String toolDirectory, String descriptionDirectory) {
        // Directory where runnable files are stored
        this.toolDirectory = toolDirectory;
        
        // Directory where application descriptions are stored
        this.descriptionDirectory = descriptionDirectory;
    }

    @Override
    protected void execute() throws Exception {
      
        // Get parameter values from user's input (order is significant)
        inputParameters = new LinkedList<String>(inputMessage.getParameters());
        
        // Get the ACD description
        acdDescription = getACD();
        
        // Get parameter descriptions from ACD file
        LinkedList<ACDParameter> acdParameters = acdDescription.getParameters();
        
        // Prepare variable map with values filled in by client
        // ACD format allows $(varname) declarations where the value
        // of varname might be not known prior to user input.
        // TODO check if we actually need this anywhere
        HashMap<String, String> varMap = new HashMap<String, String>();
        Integer index = 0;
        for (ParameterDescription param : analysis.getParameters()) {
            varMap.put(param.getName(), inputParameters.get(index));
            index++;
        }
        
        // Map description parameters to ACD parameters
        HashMap<String, ACDParameter> analysisToACD = new HashMap<String, ACDParameter>();
        ACDParameter mapTo;
        for (ParameterDescription param : analysis.getParameters()) {
            mapTo = null;
            
            // Find the appropriate ACD parameter
            for (ACDParameter acdParam : acdParameters) {
                if (acdParam.getName().equals(param.getName())) {
                    mapTo = acdParam;
                    break;
                }
            }
            analysisToACD.put(param.getName(), mapTo);
        }
        
        // Generate qualifiers and validate them
        index = 0;
        qualifiers = new LinkedList<EmbossQualifier>();
        for (ParameterDescription param : analysis.getParameters()) {
            EmbossQualifier qualifier = new EmbossQualifier(analysisToACD.get(param.getName()),
                                                            inputParameters.get(index));
            if (qualifier.isValid()) {
                qualifiers.add(qualifier);
            } else {
                // TODO Inform the user
                updateState(JobState.FAILED_USER_ERROR, "Incorrect field value: " + param.getName(), false);
                return;
            }
            index++;
        }
        
        // Processing...
        String cmd = commandLine();
        Process p = Runtime.getRuntime().exec(cmd, null, jobWorkDir);
        p.waitFor();
        
        // Some information from error stream
        byte[] b = new byte[150];
        p.getErrorStream().read(b);
        String outputString = new String(b);
        logger.info("Emboss application has finished with exit code " + p.exitValue() + 
                    " and this message: " + "\"" + outputString + "\".");
        
        // If the exit code is non-zero, the application was not successful
        if (p.exitValue() != 0) {
            logger.debug("There was an error while running emboss \"" +
                         analysis.getName() + "\" application.");
            updateState(JobState.FAILED, "EMBOSS application failed.", false);
        } 
        
        // This is what we should produce as output
        ResultMessage outputMessage = this.outputMessage;
        
        // This is where results are returned 
        ResultCallback resultHandler = this.resultHandler;
        
        outputMessage.setState(JobState.COMPLETED);
        resultHandler.sendResultMessage(inputMessage, outputMessage);
    }
    
    /**
     * Parse the appropriate ACD description file.
     * 
     * @return ACD description object.
     */
    protected ACDDescription getACD() {
        String appName = analysis.getName();
        return new ACDDescription(new File(descriptionDirectory, appName + ".acd"));
    }

    @Override
    protected void preExecute() throws Exception {
        super.preExecute();
    }
    
    /**
     * Return a command line, including executable and
     * parameters, that will be executed.
     * 
     * @return
     */
    private String commandLine() {
        
        // Form the parameters (including the executable)
        LinkedList<String> params = new LinkedList<String>();
        params.add(new File(toolDirectory, analysis.getName()).getAbsolutePath());
        
        // Parameters
        for (EmbossQualifier qualifier : qualifiers) {
            params.add(qualifier.toString());
        }
        
        // Inputs
        for (String name : inputMessage.payloadNames()) {
            params.add("-" + name + " " + name);
        }
        
        // Outputs
        // TODO: add outputs according to ACD
        params.add("-outfile " +  "outfile");
               
        String cmd = "";
        for (String string : params) {
            cmd += string + " ";
        }
        
        return cmd;
    }
    
    /**
     * Defines a qualifier which will later be passed as
     * a command line argument e.g. <pre>appname -qualname qualvalue</pre>
     * 
     * @author naktinis
     *
     */
    class EmbossQualifier {
        
        ACDParameter acdParameter;
        String value;
        
        public EmbossQualifier(ACDParameter acdParameter, String value) {
            this.acdParameter = acdParameter;
            this.value = value;
        }
        
        public boolean isValid() {
            return acdParameter.validate(value);
        }
        
        public String toString() {
            return "-" + acdParameter.getName() + " " + value;
        }
    }
}
