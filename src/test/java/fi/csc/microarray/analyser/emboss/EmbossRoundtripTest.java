package fi.csc.microarray.analyser.emboss;

import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.net.URL;
import java.util.HashMap;
import java.util.Random;

import javax.jms.JMSException;

import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;

import fi.csc.microarray.analyser.AnalysisJob;
import fi.csc.microarray.analyser.ResultCallback;
import fi.csc.microarray.analyser.ToolDescription;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.filebroker.FileBrokerClient;
import fi.csc.microarray.filebroker.FileBrokerClientMock;
import fi.csc.microarray.filebroker.FileBrokerClient.FileBrokerArea;
import fi.csc.microarray.messaging.JobState;
import fi.csc.microarray.messaging.message.ChipsterMessage;
import fi.csc.microarray.messaging.message.JobMessage;
import fi.csc.microarray.messaging.message.ResultMessage;

public class EmbossRoundtripTest {
    
    private static String path = "src/test/resources/";
    
    private boolean isResultOK = false; 

    @BeforeTest
    protected void setUp() throws Exception {
    	DirectoryLayout.uninitialise();
        DirectoryLayout.initialiseSimpleLayout();
    }

    /**
     * Test EMBOSS analysis processing roundtrip: generate description on the compute service
     * side, create new job from the description on the client side and process results
     * on the service side. Both sides are simulated, i.e., client or compute service
     * is not actually started.
     * @throws Exception 
     */
    @Test
    public void testRoundtripValidation() throws Exception {       
        // Client processes SADL and generates new analysis 
        // job according to user's input
        JobMessage jobMessage = new JobMessage();
        jobMessage.setJobId("water-" + new Random().nextInt(1000));
        jobMessage.addParameter("110.0");  // gapopen (error: greater than maximum)
        jobMessage.addParameter("10.0");   // gapextend
        jobMessage.addParameter("Y");      // brief
        
        // Process the job at compute server side
        executeJob("water.acd", jobMessage);
        
        // Check if execution failed
        Assert.assertFalse(isResultOK);
    }
    
    @Test
    public void testRoundtripExecution() throws Exception {
        // Client processes SADL and generates new analysis 
        // job according to user's input
        JobMessage jobMessage = new JobMessage();
        jobMessage.setJobId("water-" + new Random().nextInt(1000));
        jobMessage.addParameter("100.0");  // gapopen
        jobMessage.addParameter("10.0");   // gapextend
        jobMessage.addParameter("Y");      // brief
        
        // User uploads two files for input
        InputStream firstInput = new FileInputStream(path + "sequences/human_adh6.fasta");
        InputStream secondInput = new FileInputStream(path + "sequences/funghi_adh6.fasta");
        URL firstUrl = resultCallback.getFileBrokerClient().addFile(FileBrokerArea.CACHE, firstInput, -1, null);
        URL secondUrl = resultCallback.getFileBrokerClient().addFile(FileBrokerArea.CACHE, secondInput, -1, null);
        jobMessage.addPayload("asequence", firstUrl);
        jobMessage.addPayload("bsequence", secondUrl);
        
        // Process the job at compute server side
        executeJob("water.acd", jobMessage);
        
        // Check if result is ok
        Assert.assertTrue(isResultOK);
    }
    
    /**
     * Simulate compute server.
     * 
     * @param sadl
     * @param jobMessage
     */
    private void executeJob(String acdFileName, JobMessage jobMessage) throws Exception {
        // Imitate a configuration file
        HashMap<String, String> params = new HashMap<String, String>();
        params.put("externalToolPath", "/opt/chipster/tools");
        params.put("toolPath", "/EMBOSS-6.2.0/emboss");
        params.put("descriptionPath", "/EMBOSS-6.2.0/emboss/acd");
        
        // Create a job using a handler
        EmbossAnalysisHandler analysisHandler = new EmbossAnalysisHandler(params);
        ToolDescription description = analysisHandler.handle(null, acdFileName, new HashMap<String, String>()); // module should not be null
        AnalysisJob analysisJob = analysisHandler.createAnalysisJob(jobMessage,
                                                                    description, resultCallback);
        analysisJob.run();
    }
    
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

        public void sendResultMessage(ChipsterMessage inputMessage, ResultMessage resultMessage) {
            if (resultMessage.getState() == JobState.COMPLETED) {
                // This is a bit ugly and might cause trouble when tests are run in parallel?
                isResultOK = true;
            }
        }

        public boolean shouldSweepWorkDir() {
            return true;
        }
        
    };

    public static void main(String[] args) throws Exception {
        EmbossRoundtripTest test = new EmbossRoundtripTest();
        test.setUp();
        test.testRoundtripValidation();
        test.testRoundtripExecution();
        System.exit(0);
    }
    
}
