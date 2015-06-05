package fi.csc.microarray.analyser.emboss;

import java.io.File;
import java.util.HashMap;
import java.util.LinkedList;

import org.apache.log4j.Logger;

import fi.csc.microarray.analyser.JobCancelledException;
import fi.csc.microarray.analyser.ToolDescription.ParameterDescription;
import fi.csc.microarray.analyser.shell.ShellAnalysisJobBase;
import fi.csc.microarray.analyser.shell.ShellAnalysisJob.ShellParameterSecurityPolicy;
import fi.csc.microarray.messaging.JobState;
import fi.csc.microarray.messaging.message.JobMessage.ParameterValidityException;
import fi.csc.microarray.util.Exceptions;

/**
 * Runs EMBOSS applications.
 * 
 * @author naktinis, akallio, hupponen
 *
 */
public class EmbossAnalysisJob extends ShellAnalysisJobBase {
    
	public static class EmbossParameterSecurityPolicy extends ShellParameterSecurityPolicy {
		
		public boolean isValueValid(String value, ParameterDescription parameterDescription) {
			
			// Check parameter size (DOS protection) and content (shell code injection).
			// We don't check the parameter for ACD expression
			// vulnerabilities, because due to recursion limited
			// parsing and very limited syntax there should be none.
			return super.isValueValid(value, parameterDescription);
		}

	}
	
	public static EmbossParameterSecurityPolicy EMBOSS_PARAMETER_SECURITY_POLICY = new EmbossParameterSecurityPolicy();

    private String toolDirectory;
    private String descriptionDirectory;

    static final Logger logger = Logger.getLogger(EmbossAnalysisJob.class);
    
    // Output formats specified by user
    private HashMap<String, String> outputFormats = new HashMap<String, String>();
    private LinkedList<EmbossQualifier> qualifiers;
    private ACDDescription acdDescription;
    
    
    public EmbossAnalysisJob(String toolDirectory, String descriptionDirectory) {
        // Directory where runnable files are stored
        this.toolDirectory = toolDirectory;
        
        // Directory where application descriptions are stored
        this.descriptionDirectory = descriptionDirectory;
    }
    
    
    @Override
    protected void preExecute() throws JobCancelledException {
        super.preExecute();
        
        // Get parameter values from user's input (order is significant)
        LinkedList<String> inputParameters;
		try {
			inputParameters = new LinkedList<String>(inputMessage.getParameters(EMBOSS_PARAMETER_SECURITY_POLICY, analysis));
			
		} catch (ParameterValidityException e) {
			outputMessage.setErrorMessage(e.getMessage()); // always has a message
			outputMessage.setOutputText(Exceptions.getStackTrace(e));
			updateState(JobState.FAILED_USER_ERROR, "");
			return;
		}
        
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

        // store the command for execute()
        command = commandLine();
    }


    /**
     * Parse the appropriate ACD description file.
     * 
     * @return ACD description object.
     */
    private ACDDescription getACD() {
        String appName = analysis.getID();
        return new ACDDescription(new File(descriptionDirectory, appName));
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
                params.add(SHELL_STRING_SEPARATOR + qualifier.getValue() + SHELL_STRING_SEPARATOR);
            }
        }
        
        // Inputs
        for (String name : inputMessage.getKeys()) {
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
        
        // Turn off prompts
        params.add("-auto");
        
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
