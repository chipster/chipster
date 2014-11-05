package fi.csc.microarray.filebroker;

public class MockJettyFileServer extends JettyFileServer {

	public MockJettyFileServer() {
		super(null, null, null);
	}
	
	@Override
	public void start(String resourceBase, int port, String protocol) throws Exception {
		// empty
	}
	
	@Override	
	public boolean isRunning() {
		return true;
	}

}
