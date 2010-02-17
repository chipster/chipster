package fi.csc.microarray.analyser.emboss;

import java.io.File;

import org.testng.Assert;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.Test;

import fi.csc.microarray.analyser.AnalysisDescription;
import fi.csc.microarray.analyser.AnalysisDescriptionGenerator;
import fi.csc.microarray.analyser.AnalysisJob;
import fi.csc.microarray.analyser.ResultCallback;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.description.SADLDescription.Parameter;
import fi.csc.microarray.filebroker.FileBrokerClient;
import fi.csc.microarray.messaging.message.JobMessage;
import fi.csc.microarray.messaging.message.NamiMessage;
import fi.csc.microarray.messaging.message.ResultMessage;

public class EmbossRoundtripTest {
	
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
	 */
	@Test
	public void testRoundtrip() {
    
		// generate SADL
		ACDDescription acd = ACDToSADLTest.getTestDescription();
        ACDToSADL converter = new ACDToSADL(acd);
        SADLDescription sadl = converter.convert();

        // imitate client that processes SADL and generates new 
        // analysis job according to user's input
        JobMessage jobMessage = new JobMessage();
        
        for (Parameter parameter : sadl.parameters()) {
        	String paramValue = parameter.getDefaultValue(); // use defaults for all parameters
        	jobMessage.addParameter(paramValue);
        }
        // TODO we should also fill in payloads (inputs) to jobMessage
        
        // process the job at compute server side
        AnalysisDescription analysisDescription = new AnalysisDescriptionGenerator().generate(sadl, null); // FIXME we get NPE, because type is null in parameter "brief"
        EmbossAnalysisJob analysisJob = new EmbossAnalysisJob();
        analysisJob.construct(jobMessage, analysisDescription, resultCallback);
        
        // check that result is ok
        Assert.assertTrue(isResultOK);
    
    }
    
    public static void main(String[] args) throws Exception {
    	EmbossRoundtripTest test = new EmbossRoundtripTest();
		test.setUp();
		test.testRoundtrip();
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
			// TODO verify resultMessage here
			
			isResultOK = true; // this is a bit ugly and might cause trouble when tests are run in parallel?
		}

		public boolean shouldSweepWorkDir() {
			return false;
		}
		
	};

}
