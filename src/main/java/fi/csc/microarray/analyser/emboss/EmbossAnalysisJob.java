package fi.csc.microarray.analyser.emboss;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.concurrent.CountDownLatch;

import org.apache.log4j.Logger;

import fi.csc.microarray.analyser.JobCancelledException;
import fi.csc.microarray.analyser.OnDiskAnalysisJobBase;
import fi.csc.microarray.analyser.ResultCallback;
import fi.csc.microarray.analyser.AnalysisDescription.ParameterDescription;
import fi.csc.microarray.messaging.JobState;
import fi.csc.microarray.messaging.message.ResultMessage;
import fi.csc.microarray.util.Strings;

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
    
    // Output formats specified by user
    HashMap<String, String> outputFormats = new HashMap<String, String>();
    
    // Operating system process
    private Process process = null;
    
    // Latch for cancelling or finishing a job
    private CountDownLatch latch = new CountDownLatch(1);
    
    public EmbossAnalysisJob(String toolDirectory, String descriptionDirectory) {
        // Directory where runnable files are stored
        this.toolDirectory = toolDirectory;
        
        // Directory where application descriptions are stored
        this.descriptionDirectory = descriptionDirectory;
    }
    
    /**
     * User decided to cancel this job.
     */
    @Override
    protected void cancelRequested() {
        latch.countDown();
    }
    
    /**
     * Run EMBOSS job as an operating system process.
     */
    @Override
    protected void execute() throws JobCancelledException {
        
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
            // Handle non-acd parameters separately
            if (param.getName().startsWith(ACDToSADL.OUTPUT_TYPE_PREFIX)) {
                String value = !(inputParameters.get(index).equals(ACDParameter.UNDEFINED))?
                        inputParameters.get(index) : null;
                outputFormats.put(
                        param.getName().substring(ACDToSADL.OUTPUT_TYPE_PREFIX.length()),
                        value);
                continue;
            }
            
            // Handle normal parameters
            EmbossQualifier qualifier = new EmbossQualifier(analysisToACD.get(param.getName()),
                                                            inputParameters.get(index));
            EmbossQualifier.ValidityCheck check = qualifier.validate();
            if (check.passed()) {
                qualifiers.add(qualifier);
            } else {
                // Inform the user
                outputMessage.setErrorMessage(check.getMessage());
                updateState(JobState.FAILED_USER_ERROR, "Incorrect field value: " + param.getName());
                logger.debug(check.getMessage());
                return;
            }
            index++;
        }
        
        // Processing...
        try {
            String[] cmd = commandLine();
            
            // Log that we are about to run this
            logger.info("Running Emboss application " + cmd[0]);
            logger.info("Parameters: " + Strings.delimit(Arrays.asList(cmd), " "));
            
            // Start a process on the operating system
            process = Runtime.getRuntime().exec(cmd, null, jobWorkDir);
            
            // Start a new thread to listen to OS process status
            new ProcessWaiter(process, latch).start();
            
            // Job finished successfully or was cancelled
            latch.await();
            
            // Some information from error stream
            BufferedReader bufferedReader = new BufferedReader(
                    new InputStreamReader(process.getErrorStream()));
            StringBuilder stringBuilder = new StringBuilder();
            String line = null;
            while ((line = bufferedReader.readLine()) != null) {
                stringBuilder.append(line + "\n");
            }
            bufferedReader.close();
            String outputString = stringBuilder.toString();          
            logger.info("Emboss application has finished with exit code " +
                        process.exitValue() + " and this message: " +
                        "\"" + outputString + "\".");
            
            // If the exit code is non-zero, the application was not successful
            if (process.exitValue() != 0) {
                logger.debug("There was an error while running emboss \"" +
                             analysis.getDisplayName() + "\" application.");
                outputMessage.setErrorMessage(outputString);
                updateState(JobState.FAILED, "EMBOSS application failed.");
            } 
            
            // This is what we should produce as output
            ResultMessage outputMessage = this.outputMessage;
            
            // This is where results are returned 
            ResultCallback resultHandler = this.resultHandler;
            
            outputMessage.setState(JobState.RUNNING);
            resultHandler.sendResultMessage(inputMessage, outputMessage);
        } catch (IOException e) {
            // Program was not found
            logger.debug("There was an error while running emboss \"" +
                    analysis.getDisplayName() + "\" application.");
            outputMessage.setErrorMessage("Program " + acdDescription.getName() +
                    " couldn't be started.");
            updateState(JobState.FAILED, "EMBOSS application failed.");
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }
    
    /**
     * Parse the appropriate ACD description file.
     * 
     * @return ACD description object.
     */
    protected ACDDescription getACD() {
        String appName = analysis.getID();
        return new ACDDescription(new File(descriptionDirectory, appName));
    }
    
    /**
     * Destroy operating system process if it is still
     * running.
     */
    @Override
    protected void cleanUp() {
        super.cleanUp();
        process.destroy();
    }

    @Override
    protected void preExecute() throws JobCancelledException {
        super.preExecute();
    }
    
    /**
     * Return a command line that will be executed,
     * including executable and parameters.
     * 
     * @return
     */
    private String[] commandLine() {
        
        // Form the parameters (including the executable)
        LinkedList<String> params = new LinkedList<String>();

        params.add(new File(toolDirectory, analysis.getDisplayName()).getAbsolutePath());
        
        // Parameters
        for (EmbossQualifier qualifier : qualifiers) {
            if (!qualifier.getValue().equals("")) {
                params.add("-" + qualifier.getName());
                params.add(qualifier.getValue());
            }
        }
        
        // Inputs
        for (String name : inputMessage.payloadNames()) {
            params.add("-" + name);
            params.add(name);
        }
        
        // Simple outputs
        for (ACDParameter param : acdDescription.getOutputParameters()) {
            // User might have changed output type
            String format = outputFormats.get(param.getName());
            if (format != null) {
                // Add additional qualifier for defining format
                if (param.getType().equals("align")) {
                    params.add("-aformat");
                } else {
                    params.add("-osformat");
                }
                params.add(format);
            }
            
            // Add qualifier
            params.add("-" + param.getName());
            params.add(param.getOutputFilename(true));
        }
        
        // Graphics outputs
        for (ACDParameter param : acdDescription.getGraphicsParameters()) {
            // Emboss automatically adds extensions for graphics files
            params.add("-" + param.getName());
            params.add("png");
            params.add("-goutfile");
            params.add(param.getOutputFilename(false));
        }
        
        String[] cmd = new String[0];
        return params.toArray(cmd);
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
        
        public String getName() {
            return acdParameter.getName();
        }
        
        public String getValue() {
            return acdParameter.normalize(value);
        }
        
        public ValidityCheck validate() {
            // Check if this value is ok for this parameter
            if (acdParameter.isRequired() && value.equals("")) {
                return new ValidityCheck(false, "Parameter \"" + acdParameter.getName() +
                                                "\" can not be empty.");
            } else if (!value.equals("") && !acdParameter.validate(value)) {
                return new ValidityCheck(false, "Incorrect value \"" + getValue() +
                                                "\" for \"" + acdParameter.getName() + "\" parameter");
            }
            return new ValidityCheck(true, null);
        }
    }
}
