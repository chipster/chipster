package fi.csc.microarray.messaging;

import java.io.IOException;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URL;

import org.testng.annotations.Test;

public class HttpTest {

	
	@Test(groups = {"smoke"} )
	public void URLCheck() throws MalformedURLException, IOException {
		HttpURLConnection con;
		URL url = new URL("http://ocicat.csc.fi:8080/fileserver/repository/" + "004341ec-af2d-429d-86ed-a74e4c845a75_normalized.tsv");
		System.out.println(url);
		con = (HttpURLConnection) url.openConnection();
		System.out.println(con.getResponseCode());
		System.out.println(con.getContentLength());
		System.out.println(Integer.MAX_VALUE);
	
	}
}
