package fi.csc.microarray.analyser.emboss;

import java.io.File;
import java.io.FileInputStream;

import org.testng.Assert;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.Test;

import fi.csc.microarray.analyser.AnalysisDescription;
import fi.csc.microarray.analyser.AnalysisJob;
import fi.csc.microarray.analyser.ResultCallback;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.fs.FSDataManager;
import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.filebroker.FileBrokerClient;
import fi.csc.microarray.messaging.JobState;
import fi.csc.microarray.messaging.message.JobMessage;
import fi.csc.microarray.messaging.message.NamiMessage;
import fi.csc.microarray.messaging.message.ResultMessage;

public class EmbossRoundtripTest {
    
    private static String path = "src/test/resources/";
    
    private boolean isResultOK = false; 

    @BeforeSuite
    protected void setUp() throws Exception {
        DirectoryLayout.initialiseUnitTestLayout();            
    }

    /**
     * Test EMBOSS analysis processing roundtrip: generate description on the compute service
     * side, create new job from the description on the client side and process results
     * on the service side. Both sides are simulated, i.e., client or compute service
     * is not actually started.
     * @throws Exception 
     */
    @Test
    public void testRoundtrip() throws Exception {
        // Client processes SADL and generates new analysis 
        // job according to user's input
        JobMessage jobMessage = new JobMessage();
        jobMessage.addParameter("110.0");  // gapopen (error: greater than maximum)
        jobMessage.addParameter("10.0");   // gapextend
        jobMessage.addParameter("Y");      // brief
        // TODO test input files
        //DataBean bean = beanFromFile("sequences/human_adh6.fasta"); 
        //jobMessage.addPayload("asequence", bean.getUrl());
        
        // Process the job at compute server side
        executeJob("water.acd", jobMessage);
        
        // Check that result is ok
        Assert.assertTrue(!isResultOK);
    }
    
    public static void main(String[] args) throws Exception {
        EmbossRoundtripTest test = new EmbossRoundtripTest();
        test.setUp();
        test.testRoundtrip();
    }
    
    private DataBean beanFromFile(String filename) {
        FileInputStream input;
        try {
            input = new FileInputStream(new File(path + filename));
            DataBean bean = new FSDataManager().createDataBean(filename, input);
            return bean;
        } catch (Exception e) {
            e.printStackTrace();
        }
        return null;
    }
    
    /**
     * Create SADL description from given EMBOSS application name.
     * 
     * @param appName
     * @return
     */
    private SADLDescription createSADL(String appName) {
        ACDDescription acd = ACDToSADLTest.getTestDescription("water");
        ACDToSADL converter = new ACDToSADL(acd);
        return converter.convert();
    }
    
    /**
     * Simulate compute server.
     * 
     * @param sadl
     * @param jobMessage
     */
    private void executeJob(String acdFileName, JobMessage jobMessage) throws Exception {
        EmbossAnalysisHandler analysisHandler = new EmbossAnalysisHandler();
        AnalysisDescription description = analysisHandler.handle(acdFileName);
        AnalysisJob analysisJob = analysisHandler.createAnalysisJob(jobMessage,
                                                                    description, resultCallback);
        analysisJob.run();
    }
    
    private ResultCallback resultCallback = new ResultCallback() {

        public FileBrokerClient getFileBrokerClient() {
            return null;
        }

        public File getWorkDir() {
            return null;
        }

        public void removeRunningJob(AnalysisJob job) {
        }

        public void sendResultMessage(NamiMessage inputMessage, ResultMessage resultMessage) {          
            if (resultMessage.getState() == JobState.COMPLETED) {
                // This is a bit ugly and might cause trouble when tests are run in parallel?
                isResultOK = true;
            }
        }

        public boolean shouldSweepWorkDir() {
            return false;
        }
        
    };

}
