package fi.csc.microarray.analyser.shell;

import java.util.LinkedList;

import fi.csc.microarray.analyser.JobCancelledException;
import fi.csc.microarray.analyser.AnalysisDescription.InputDescription;
import fi.csc.microarray.analyser.AnalysisDescription.OutputDescription;
import fi.csc.microarray.analyser.AnalysisDescription.ParameterDescription;

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
        Boolean inputLast = analysis.getConfigParameters().get("input") != null &&
                analysis.getConfigParameters().get("input").toLowerCase().equals("last");
        
        // Additional arguments
        String arguments = analysis.getConfigParameters().get("arguments");
        String[] extraArguments = (arguments != null && !arguments.equals("")) ?
                arguments.split(",") : new String[] {};
 
        String outputParameter = null;
        if (!useStdout) {    
            // If program creates a normal file, we need an output parameter
            outputParameter = analysis.getConfigParameters().get("output");
        }

        LinkedList<String> inputParameters;
    	
    	// Get parameter values from user's input (order is significant)
        inputParameters = new LinkedList<String>(inputMessage.getParameters());
                
        // Generate the command to be executed
        LinkedList<String> commandParts = new LinkedList<String>();
        commandParts.add(executablePath);
        
        // Prepend arguments defined in the configuration file
        for (String arg : extraArguments) {
            commandParts.add(arg);
        }
        
        // Parameters
        int index = 0;
        for (ParameterDescription parameter : analysis.getParameters()) {
            String value = inputParameters.get(index);
            if (!value.equals("")) {
                commandParts.add("-" + parameter.getName());
                commandParts.add(value);
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
        for (InputDescription input : analysis.getInputFiles()) {
        	if (!inputLast) {
                // Input is a named parameter
                commandParts.add("-" + input.getFileName());
            }
            commandParts.add(input.getFileName());
        }
        
        // store the command for execute()
        command = commandParts.toArray(new String[] {});
    }
}
