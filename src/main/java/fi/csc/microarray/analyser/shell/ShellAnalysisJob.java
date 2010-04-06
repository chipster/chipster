package fi.csc.microarray.analyser.shell;

import fi.csc.microarray.analyser.JobCancelledException;
import fi.csc.microarray.analyser.OnDiskAnalysisJobBase;

/**
 * Job that is run as a generic shell command.
 * 
 * @author naktinis
 *
 */
public class ShellAnalysisJob extends OnDiskAnalysisJobBase {
    
    private String toolDirectory;
    private String descriptionDirectory;
    
    public ShellAnalysisJob(String toolDirectory, String descriptionDirectory) {
        // Directory where runnable files are stored
        this.toolDirectory = toolDirectory;
        
        // Directory where application descriptions are stored
        this.descriptionDirectory = descriptionDirectory;
        System.out.println(toolDirectory);
    }

    @Override
    protected void cancelRequested() { }

    @Override
    protected void execute() throws JobCancelledException {
        // TODO Auto-generated method stub
        
    }

}
