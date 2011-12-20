package fi.csc.microarray.analyser.shell;

import java.util.LinkedList;

import fi.csc.microarray.analyser.JobCancelledException;
import fi.csc.microarray.analyser.AnalysisDescription.InputDescription;
import fi.csc.microarray.analyser.AnalysisDescription.OutputDescription;
import fi.csc.microarray.analyser.AnalysisDescription.ParameterDescription;
import fi.csc.microarray.messaging.JobState;
import fi.csc.microarray.messaging.message.JobMessage.ParameterSecurityPolicy;
import fi.csc.microarray.messaging.message.JobMessage.ParameterValidityException;
import fi.csc.microarray.util.Exceptions;

/**
 * Job that is run as a generic shell command.
 * <p>
 * Parameters that could be given in the tool configuration:
 *  <ul>
 *  <li> output - name of the parameter for passing output files
 *  <li> stdout - if equals to "yes", the output is read from stdout
 *  <li> input - if equals to "last", then input is given without parameter
 *               name and as the last parameter
 *  <li> arguments - comma-separated list of arguments that will get passed
 *               as the first arguments for the command.   
 *  </ul>
 * 
 * @author naktinis, hupponen
 *
 */
public class ShellAnalysisJob extends ShellAnalysisJobBase {
    
	public static class ShellParameterSecurityPolicy implements ParameterSecurityPolicy {
		
		private static final int MAX_VALUE_LENGTH = 1000;
		public static String COMMAND_LINE_SAFE_VALUE_PATTERN = "[\\w+\\-_:\\.,*()= ]*"; // Only word characters and some special symbols are allowed
	
		public boolean isValueValid(String value, ParameterDescription parameterDescription) {
			
			// Check parameter size (DOS protection)
			if (value.length() > MAX_VALUE_LENGTH) {
				return false;
			}
			
			// Check content for string termination (shell code injection)
			if (value.contains(SHELL_STRING_SEPARATOR)) {
				return false;
			}

			// Check that content matches the pattern (shell code injection)
			return value.matches(COMMAND_LINE_SAFE_VALUE_PATTERN);
		}

	}
	
	public static ShellParameterSecurityPolicy SHELL_PARAMETER_SECURITY_POLICY = new ShellParameterSecurityPolicy();

	private static final String USE_ONLY_PARAMETER_VALUE_TOKEN = "-value-only";
	private static final String NO_PARAMETER_VALUE_TOKEN = "NO-VALUE";
	
    /**
     * Construct the command line.
     * 
     */
    @Override
    protected void preExecute() throws JobCancelledException {
    	cancelCheck();
    	super.preExecute();
    	
        // Path to executable file
        String executablePath = analysis.getCommand();
        
        // Output method
        useStdout = analysis.getConfigParameters().get("stdout") != null &&
                analysis.getConfigParameters().get("stdout").toLowerCase().equals("yes");
        
        // Input parameter
        boolean inputLast = analysis.getConfigParameters().get("input") != null &&
                analysis.getConfigParameters().get("input").toLowerCase().equals("last");
        
        boolean noInput = analysis.getConfigParameters().get("input") != null &&
        analysis.getConfigParameters().get("input").toLowerCase().equals("none");
        
        boolean inputsAsPlainArguments = analysis.getConfigParameters().get("inputs-as-plain-arguments") != null &&
        analysis.getConfigParameters().get("inputs-as-plain-arguments").toLowerCase().equals("true");
        
        
        // Additional arguments
        String arguments = analysis.getConfigParameters().get("arguments");
        String[] extraArguments = (arguments != null && !arguments.equals("")) ?
                arguments.trim().split(" ") : new String[] {};
 
        String outputParameter = null;
        if (!useStdout) {    
            // If program creates a normal file, we need an output parameter
            outputParameter = analysis.getConfigParameters().get("output");
        }

        LinkedList<String> inputParameters;
    	
    	// Get parameter values from user's input (order is significant)
        try {
			inputParameters = new LinkedList<String>(inputMessage.getParameters(SHELL_PARAMETER_SECURITY_POLICY, analysis));
		} catch (ParameterValidityException e) {
			outputMessage.setErrorMessage(e.getMessage()); // always has a message
			outputMessage.setOutputText(Exceptions.getStackTrace(e));
			updateState(JobState.FAILED_USER_ERROR, "");
			return;
		}
                
        // Generate the command to be executed
        LinkedList<String> commandParts = new LinkedList<String>();
        commandParts.add(executablePath);

        // if plain inputs, add them right after the executable
        if (inputsAsPlainArguments) {
        	for (InputDescription input : analysis.getInputFiles()) {
        		commandParts.add(input.getFileName());
        	}
        }
        
        // Prepend arguments defined in the configuration file
        for (String arg : extraArguments) {
            commandParts.add(arg);
        }
        
        // Parameters
        int index = 0;
        for (ParameterDescription parameter : analysis.getParameters()) {
            String value = inputParameters.get(index);

            // value only parameters
            if (parameter.getName().endsWith(USE_ONLY_PARAMETER_VALUE_TOKEN)) {
            	
            	// no value parameter, don't add anything
            	if (!value.equals(NO_PARAMETER_VALUE_TOKEN)) {
            		commandParts.add(value);
//            		commandParts.add(SHELL_STRING_SEPARATOR + value + SHELL_STRING_SEPARATOR);
            	}
            	
            }
            
            
            // normal parameters
            else if (!value.equals("")) {
                commandParts.add("-" + parameter.getName());
                commandParts.add(value);
//              commandParts.add(SHELL_STRING_SEPARATOR + value + SHELL_STRING_SEPARATOR);
            }
            index++;
        }

        // Outputs to a file (currently we only support a single output)
        if (outputParameter != null) {
            OutputDescription output = analysis.getOutputFiles().get(0);
            commandParts.add("-" + outputParameter);
            commandParts.add(output.getFileName().getID());
        }
        
        // Inputs
        if (!inputsAsPlainArguments && !noInput) {
        	for (InputDescription input : analysis.getInputFiles()) {
        		if (!inputLast) {
        			// Input is a named parameter
        			commandParts.add("-" + input.getFileName());
        		}
        		commandParts.add(input.getFileName());
        	}
        }
        	
        // store the command for execute()
        command = commandParts.toArray(new String[] {});
    }
}
