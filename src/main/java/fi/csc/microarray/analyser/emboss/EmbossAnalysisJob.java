package fi.csc.microarray.analyser.emboss;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
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
        HashMap<String, String> varMap = new HashMap<String, String>();
        Integer index = 0;
        for (ParameterDescription param : analysis.getParameters()) {
            varMap.put(param.getName(), inputParameters.get(index));
            index++;
        }
        acdDescription.updateVariables(varMap);
        
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
            EmbossQualifier.ValidityCheck check = qualifier.validate();
            if (check.passed()) {
                qualifiers.add(qualifier);
            } else {
                // Inform the user
                outputMessage.setErrorMessage(check.getMessage());
                updateState(JobState.FAILED_USER_ERROR, "Incorrect field value: " + param.getName(), false);
                logger.debug(check.getMessage());
                return;
            }
            index++;
        }
        
        // Processing...
        String cmd = commandLine();
        Process p = Runtime.getRuntime().exec(cmd, null, jobWorkDir);
        p.waitFor();
        
        // Some information from error stream
        BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(p.getErrorStream()));
        StringBuilder stringBuilder = new StringBuilder();
        String line = null;
        while ((line = bufferedReader.readLine()) != null) {
            stringBuilder.append(line + "\n");
        }
        bufferedReader.close();
        String outputString = stringBuilder.toString();
        
        logger.info("Running Emboss application " + cmd);
        logger.info("Emboss application has finished with exit code " + p.exitValue() + 
                    " and this message: " + "\"" + outputString + "\".");
        
        // If the exit code is non-zero, the application was not successful
        if (p.exitValue() != 0) {
            logger.debug("There was an error while running emboss \"" +
                         analysis.getName() + "\" application.");
            outputMessage.setErrorMessage(outputString);
            updateState(JobState.FAILED, "EMBOSS application failed.", false);
        } 
        
        // This is what we should produce as output
        ResultMessage outputMessage = this.outputMessage;
        
        // This is where results are returned 
        ResultCallback resultHandler = this.resultHandler;
        
        outputMessage.setState(JobState.RUNNING);
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
        
        // Simple outputs
        for (ACDParameter param : acdDescription.getOutputParameters()) {
            params.add("-" + param.getName() + " " + param.getOutputFilename(true));
        }
        
        // Graphics outputs
        for (ACDParameter param : acdDescription.getGraphicsParameters()) {
            // Emboss automatically adds extensions for graphics files
            params.add("-" + param.getName() + " png");
            params.add("-goutfile " + param.getOutputFilename(false));
        }
                       
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
        
        /**
         * Represents validity check. Encapsulates check
         * result and the error message if any.
         */
        class ValidityCheck {
            private boolean passed;
            private String message;
            
            public ValidityCheck(boolean passed, String message) {
                this.passed = passed;
                this.message = message;
            }
            
            public boolean passed() {
                return this.passed;
            }
            
            public String getMessage() {
                return this.message;
            }
        }
        
        ACDParameter acdParameter;
        String value;
        
        public EmbossQualifier(ACDParameter acdParameter, String value) {
            this.acdParameter = acdParameter;
            this.value = value;
        }
        
        public String getValue() {
            return value;
        }
        
        public ValidityCheck validate() {
            // Check if this value is ok for this parameter
            if (acdParameter.isRequired() && value.equals("")) {
                return new ValidityCheck(false, "Parameter \"" + acdParameter.getName() +
                                                "\" can not be empty.");
            } else if (!acdParameter.validate(value)) {
                return new ValidityCheck(false, "Incorrect value \"" + getValue() +
                                                "\" for \"" + acdParameter.getName() + "\" parameter");
            }
            return new ValidityCheck(true, null);
        }
        
        public String toString() {
            // If value is empty, don't include this qualifier at all
            if (!value.equals("")) {
                return "-" + acdParameter.getName() + " " + value;
            }
            return "";
        }
    }
}
