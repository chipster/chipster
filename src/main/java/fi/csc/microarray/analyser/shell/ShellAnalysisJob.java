package fi.csc.microarray.analyser.shell;

import java.util.LinkedList;

import org.apache.log4j.Logger;

import fi.csc.microarray.analyser.AnalysisDescription;
import fi.csc.microarray.analyser.JobCancelledException;
import fi.csc.microarray.analyser.OnDiskAnalysisJobBase;
import fi.csc.microarray.analyser.ResultCallback;
import fi.csc.microarray.analyser.AnalysisDescription.ParameterDescription;
import fi.csc.microarray.analyser.emboss.EmbossAnalysisJob;
import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.description.SADLParser;
import fi.csc.microarray.description.SADLDescription.Input;
import fi.csc.microarray.description.SADLParser.ParseException;
import fi.csc.microarray.messaging.JobState;
import fi.csc.microarray.messaging.message.ResultMessage;
import fi.csc.microarray.util.Files;

/**
 * Job that is run as a generic shell command.
 * 
 * @author naktinis
 *
 */
public class ShellAnalysisJob extends OnDiskAnalysisJobBase {
   
    private AnalysisDescription description;
    private SADLDescription sadl;
    private String executablePath;
    
    LinkedList<String> inputParameters;
    
    static final Logger logger = Logger.getLogger(EmbossAnalysisJob.class);
    
    public ShellAnalysisJob(AnalysisDescription ad) {
        // Store descriptions
        this.description = ad;
        try {
            this.sadl = new SADLParser().parse(ad.getSADL());
        } catch (ParseException e) {
            e.printStackTrace();
        }
        
        // Path to executable file
        this.executablePath = ad.getCommand();
    }

    @Override
    protected void cancelRequested() { }

    @Override
    protected void execute() throws JobCancelledException {
        // Get parameter values from user's input (order is significant)
        inputParameters = new LinkedList<String>(inputMessage.getParameters());
                
        // Generate the command to be executed
        LinkedList<String> command = new LinkedList<String>();
        command.add(executablePath);
        
        // Parameters
        int index = 0;
        for (ParameterDescription parameter : description.getParameters()) {
            command.add("-" + parameter.getName());
            command.add(inputParameters.get(index));
            index++;
        }
        
        // Inputs
        for (Input input : sadl.inputs()) {
            command.add("-" + input.getName().getID());
            command.add(input.getName().getID());
        }
        
        // Outputs
        for (String output : description.getOutputFiles()) {
            command.add("-" + output);
            command.add(output);
        }
        
        String[] cmd = new String[0];
        cmd = command.toArray(cmd);
        try {
            logger.info("Running Shell application " + cmd[0]);
            
            Process p = Runtime.getRuntime().exec(cmd, null, jobWorkDir);
            p.waitFor();
            
            logger.info("Shell application has finished with exit code " + p.exitValue() + 
                        " and this message: " + "\"" +
                        Files.inputStreamToString(p.getErrorStream()) + "\".");
            
            // This is what we should produce as output
            ResultMessage outputMessage = this.outputMessage;
            
            // This is where results are returned 
            ResultCallback resultHandler = this.resultHandler;
            
            outputMessage.setState(JobState.RUNNING);
            resultHandler.sendResultMessage(inputMessage, outputMessage);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
