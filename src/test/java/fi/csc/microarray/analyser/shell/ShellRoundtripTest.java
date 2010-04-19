package fi.csc.microarray.analyser.shell;

import java.io.File;
import java.util.HashMap;
import java.util.Random;

import javax.jms.JMSException;

import org.testng.Assert;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.Test;

import fi.csc.microarray.analyser.AnalysisDescription;
import fi.csc.microarray.analyser.AnalysisJob;
import fi.csc.microarray.analyser.ResultCallback;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.filebroker.FileBrokerClient;
import fi.csc.microarray.filebroker.FileBrokerClientMock;
import fi.csc.microarray.messaging.message.ChipsterMessage;
import fi.csc.microarray.messaging.message.JobMessage;
import fi.csc.microarray.messaging.message.ResultMessage;

/**
 * Test generic command line analysis processing roundtrip: generate description on
 * the compute service side, create new job from the description on the client side
 * and process results on the service side. Both sides are simulated, i.e., client
 * or compute service is not actually started, however, the shell command is run.
 * 
 * @author naktinis
 */
public class ShellRoundtripTest {
    
    private static String path = "src/test/resources/";
    
    private boolean isResultOK = false; 

    @BeforeSuite
    protected void setUp() throws Exception {
        DirectoryLayout.initialiseClientLayout();
    }
    
    @Test
    public void testRoundtripExecution() throws Exception {
        // Client processes SADL and generates new analysis 
        // job according to user's input
        JobMessage jobMessage = new JobMessage();
        jobMessage.setJobId("import-" + new Random().nextInt(1000));
        jobMessage.addParameter("1");   // sbegin
        jobMessage.addParameter("10");  // send
        
        // Process the job at compute server side
        executeJob("import.sadl", jobMessage);
        
        // Check if result is ok
        Assert.assertTrue(isResultOK);
    }
    
    /**
     * Simulate compute server.
     * 
     * @param sadl
     * @param jobMessage
     */
    private void executeJob(String sadlFileName, JobMessage jobMessage) throws Exception {
        // Imitate a configuration file
        HashMap<String, String> params = new HashMap<String, String>();
        params.put("descriptionPath", "/opt/chipster/tools/shell");
        
        // Create a job using a handler
        ShellAnalysisHandler analysisHandler = new ShellAnalysisHandler(params);
        // Simulate parameters from sequence-module.xml
        HashMap<String, String> map = new HashMap<String, String>();
        map.put("output", "outseq");
        map.put("executable", "/opt/EMBOSS-6.2.0/emboss/seqret");
        AnalysisDescription description = analysisHandler.handle(sadlFileName, map);
        AnalysisJob analysisJob = analysisHandler.createAnalysisJob(jobMessage,
                                                                    description, resultCallback);
        analysisJob.run();
    }
    
    /**
     * Callback message that needs to be checked.
     */
    private ResultCallback resultCallback = new ResultCallback() {

        private FileBrokerClient fileBroker = null;

        public FileBrokerClient getFileBrokerClient() {
            // Create a mock file broker
            if (fileBroker == null) {
                try {
                    fileBroker = new FileBrokerClientMock();
                } catch (JMSException e) {
                    e.printStackTrace();
                }
            }
            
            return fileBroker;
        }

        public File getWorkDir() {
            // Create a temporary directory for testing purposes
            File jobDir = new File(path, "emboss-tmp");
            if (!jobDir.exists()) {
                jobDir.mkdir();
            }
            return jobDir;
        }

        public void removeRunningJob(AnalysisJob job) {
        }

        public void sendResultMessage(ChipsterMessage inputMessage,
                                      ResultMessage resultMessage) {
            // Result message has been sent
            isResultOK = true;
        }

        public boolean shouldSweepWorkDir() {
            return true;
        }
        
    };
    
    public static void main(String[] args) throws Exception {
        ShellRoundtripTest test = new ShellRoundtripTest();
        test.setUp();
        test.testRoundtripExecution();
        System.exit(0);
    }

}
