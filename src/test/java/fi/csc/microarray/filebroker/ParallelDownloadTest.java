package fi.csc.microarray.filebroker;

import java.io.IOException;
import java.io.InputStream;
import java.net.HttpURLConnection;
import java.net.URL;
import java.net.URLConnection;
import java.security.KeyManagementException;
import java.security.NoSuchAlgorithmException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import javax.jms.JMSException;

import fi.csc.microarray.util.KeyAndTrustManager;

public class ParallelDownloadTest {
	
	public static boolean readOnlyStart = true;
			
	public static void main(String args[]) throws IOException, JMSException, KeyManagementException, NoSuchAlgorithmException, InterruptedException {
		
//		System.setProperty("http.maxConnections", "10");
//		System.out.println(System.getProperty("http.maxConnections") + " connections");
		
		ExecutorService executor = Executors.newFixedThreadPool(5);
		long t = System.currentTimeMillis();
		for (int i = 0; i < 100; i++) {
			//download();
			executor.execute(new Runnable() {
				@Override
				public void run() {
					try {
						download();
					} catch (KeyManagementException | NoSuchAlgorithmException
							| IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			});
		}
		
		executor.shutdown();
		executor.awaitTermination(60, TimeUnit.SECONDS);
		System.out.println((System.currentTimeMillis() - t) + "ms");
		
		Thread.sleep(30_000);
		System.out.println("30 seconds elapsed");
	}

	public static void download() throws IOException, KeyManagementException, NoSuchAlgorithmException {
		URL url = new URL("https://vm0180.kaj.pouta.csc.fi:8080/storage/86ced148-6eba-4e38-9b45-af79d9fcc9cb");
		URLConnection connection = url.openConnection();
		KeyAndTrustManager.configureForTrustAllCertificates(connection);
		
		try (InputStream input = connection.getInputStream()) {
			long total = 0;
			while (true) {
				int bytes = input.read(new byte[64*1024]);
				total += bytes;
				if (bytes <= 0) {
					break;
				}
				if (readOnlyStart) {
					break;
				}
			}
			((HttpURLConnection)connection).disconnect();
//			System.out.println("Read " + total + " bytes");
		}
	}
}
