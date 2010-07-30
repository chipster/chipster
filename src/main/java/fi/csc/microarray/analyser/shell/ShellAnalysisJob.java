package fi.csc.microarray.analyser.shell;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.util.LinkedList;
import java.util.concurrent.CountDownLatch;

import org.apache.log4j.Logger;

import fi.csc.microarray.analyser.AnalysisDescription;
import fi.csc.microarray.analyser.JobCancelledException;
import fi.csc.microarray.analyser.OnDiskAnalysisJobBase;
import fi.csc.microarray.analyser.ResultCallback;
import fi.csc.microarray.analyser.AnalysisDescription.OutputDescription;
import fi.csc.microarray.analyser.AnalysisDescription.ParameterDescription;
import fi.csc.microarray.analyser.emboss.EmbossAnalysisJob;
import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.description.SADLParser;
import fi.csc.microarray.description.SADLDescription.Input;
import fi.csc.microarray.description.SADLParser.ParseException;
import fi.csc.microarray.messaging.JobState;
import fi.csc.microarray.messaging.message.ResultMessage;
import fi.csc.microarray.util.Files;
import fi.csc.microarray.util.Strings;

/**
 * Job that is run as a generic shell command.
 * <p>
 * Parameters that could be given in the tool configuration:
 *  <ul>
 *  <li> output - name of the parameter for passing output files
 *  <li> stdout - if equals to "yes", the output is read from stdout
 *  <li> input - if equals to "last", then input is given without parameter
 *               name and as the last parameter
 *  </ul>
 * 
 * @author naktinis
 *
 */
public class ShellAnalysisJob extends OnDiskAnalysisJobBase {
   
    private AnalysisDescription description;
    private SADLDescription sadl;
    private String executablePath;
    private Boolean useStdout;
    private Boolean inputLast;
    private String outputParameter = null;
    
    LinkedList<String> inputParameters;
    
    // Logger for this class
    static final Logger logger = Logger.getLogger(EmbossAnalysisJob.class);
    
    // Latch for cancelling or finishing a job
    private CountDownLatch latch = new CountDownLatch(1);
    
    // Operating system process
    private Process process = null;
    
    public ShellAnalysisJob(AnalysisDescription ad) {
        // Store descriptions
        this.description = ad;
        try {
            this.sadl = new SADLParser().parse(ad.getSADL());
        } catch (ParseException e) {
            e.printStackTrace();
        }
        
        // Path to executable file
        executablePath = ad.getCommand();
        
        // Output method
        useStdout = ad.getConfigParameters().get("stdout") != null &&
                ad.getConfigParameters().get("stdout").toLowerCase().equals("yes");
        
        // Input parameter
        inputLast = ad.getConfigParameters().get("input") != null &&
                ad.getConfigParameters().get("input").toLowerCase().equals("last");
        
        if (!useStdout) {    
            // If program creates a normal file, we need an output parameter
            outputParameter = ad.getConfigParameters().get("output");
        }
    }

    /**
     * User decided to cancel this job.
     */
    @Override
    protected void cancelRequested() {
        latch.countDown();
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
    protected void execute() throws JobCancelledException {
        // Get parameter values from user's input (order is significant)
        inputParameters = new LinkedList<String>(inputMessage.getParameters());
                
        // Generate the command to be executed
        LinkedList<String> command = new LinkedList<String>();
        command.add(executablePath);
        
        // Parameters
        int index = 0;
        for (ParameterDescription parameter : description.getParameters()) {
            String value = inputParameters.get(index);
            if (!value.equals("")) {
                command.add("-" + parameter.getName());
                command.add(value);
            }
            index++;
        }

        // Outputs to a file (currently we only support a single output)
        if (outputParameter != null) {
            OutputDescription output = description.getOutputFiles().get(0);
            command.add("-" + this.outputParameter);
            command.add(output.getFileName().getID());
        }
        
        // Inputs
        for (Input input : sadl.inputs()) {
            if (!inputLast) {
                // Input is a named parameter
                command.add("-" + input.getName().getID());
            }
            command.add(input.getName().getID());
        }
        
        String[] cmd = new String[0];
        cmd = command.toArray(cmd);
        try {
            logger.info("Running Shell application " + cmd[0]);
            logger.info("Parameters: " + Strings.delimit(command, " "));
            
            process = Runtime.getRuntime().exec(cmd, null, jobWorkDir);
            
            // Start a new thread to listen to OS process status
            new ProcessWaiter(process, latch).start();
            
            // Job finished successfully or was cancelled
            latch.await();
            
            String outputString = Files.inputStreamToString(process.getErrorStream());
            
            logger.info("Shell application has finished with exit code " +
                        process.exitValue() + " and this message: " +
                        "\"" + outputString + "\".");
            
            // If the exit code is non-zero, the application was not successful
            if (process.exitValue() != 0) {
                logger.debug("There was an error while running \"" +
                        analysis.getDisplayName() + "\" application.");
                outputMessage.setErrorMessage(outputString);
                updateState(JobState.FAILED, "Application failed.");
            } else if (useStdout) {
                // FIXME not tested yet
                // Take output data from stdout
                InputStream stdoutStream =
                        new BufferedInputStream(process.getInputStream());
                OutputDescription output = description.getOutputFiles().get(0);
                File outputFile = new File(jobWorkDir, output.getFileName().getID());
                FileOutputStream fileStream = new FileOutputStream(outputFile);
                
                // Write from program's stdout to output file
                byte[] buffer = new byte[4096];  
                int bytesRead;  
                while ((bytesRead = stdoutStream.read(buffer)) != -1) {  
                    fileStream.write(buffer, 0, bytesRead);  
                }  
                stdoutStream.close();
                fileStream.close();
            }
            
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
