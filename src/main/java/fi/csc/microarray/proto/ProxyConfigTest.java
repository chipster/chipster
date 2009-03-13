package fi.csc.microarray.proto;


import java.net.Proxy;
import java.net.ProxySelector;
import java.net.URI;
import java.net.URISyntaxException;

public class ProxyConfigTest {

	private static final String CHIPSTER_URL = "http://chipster.csc.fi:8080";

	public static void main(String[] args) throws Exception {
		
		System.out.println("Original proxy configuration for " + CHIPSTER_URL);
		showProxies();
		System.out.println();
		
		ProxySelector.setDefault(new NoProxySelector());
		
		System.out.println("Overwritten proxy configuration for " + CHIPSTER_URL);
		showProxies();
	}
	
	public static void showProxies() throws URISyntaxException {
		for (Proxy proxy : ProxySelector.getDefault().select(new URI(CHIPSTER_URL))) {
			System.out.println(proxy);
		}
	}
}
