package fi.csc.microarray.util.rest;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.UUID;

import org.mortbay.util.IO;
import org.testng.Assert;

public class FileServerTest {

	private String repositoryUrl;
	private String testfile;
	boolean put = true;
	boolean get = true;

	public FileServerTest(String repositoryUrl, String testfile, int threadCount, int repeatCount, boolean put, boolean get) {
		System.out.println("Rest test with repository: " + repositoryUrl + ", threads: " + threadCount + ", repeats: " + repeatCount);

		this.repositoryUrl = repositoryUrl;
		this.testfile = testfile;

		if (!put && !get) {
			this.put = true;
			this.get = true;
		} else {
			this.put = put;
			this.get = get;
		}

		for (int i = 0; i < threadCount; i++) {
			new Thread(new PutGetTest(repeatCount)).start();
		}

	}

	private boolean isSuccessfulCode(int responseCode) {
		return responseCode >= 200 && responseCode < 300; // 2xx => successful
	}

	class PutGetTest implements Runnable {

		private int repeatCount;

		public PutGetTest(int repeatCount) {
			this.repeatCount = repeatCount;
		}

		public void run() {
			String filename = "";
			boolean putOnce = false;
			if (get && !put) {
				filename = "testfile-" + UUID.randomUUID().toString();
				putOnce = true;
			}

			for (int i = 0; i < repeatCount; i++) {

				try {

					URL url = new URL(repositoryUrl + "/" + filename);
					HttpURLConnection connection = null;

					if (putOnce) {
						filename = "testfile-" + UUID.randomUUID().toString();
						url = new URL(repositoryUrl + "/" + filename);
						connection = (HttpURLConnection) url.openConnection();
						connection.setRequestMethod("PUT");
						connection.setDoOutput(true);
						connection.setChunkedStreamingMode(2048);
						OutputStream os = connection.getOutputStream();
						IO.copy(new FileInputStream(new File(testfile)), os);
						os.close();
						Assert.assertTrue(isSuccessfulCode(connection.getResponseCode()));
						connection.disconnect();
						putOnce = false;
					}

					if (put) {
						// 1. upload
						filename = "testfile-" + UUID.randomUUID().toString();
						url = new URL(repositoryUrl + "/" + filename);
						connection = (HttpURLConnection) url.openConnection();
						connection.setRequestMethod("PUT");
						connection.setDoOutput(true);
						connection.setChunkedStreamingMode(2048);
						OutputStream os = connection.getOutputStream();
						IO.copy(new FileInputStream(new File(testfile)), os);
						os.close();
						Assert.assertTrue(isSuccessfulCode(connection.getResponseCode()));
						connection.disconnect();
					}
					// Thread.sleep(1000);

					if (get) {
						// 2. download
						connection = (HttpURLConnection) url.openConnection();
						InputStream is = connection.getInputStream();
						OutputStream fos = new FileOutputStream(new File(filename));

						IO.copy(is, fos);
						is.close();
						fos.close();
						Assert.assertTrue(isSuccessfulCode(connection.getResponseCode()));
						connection.disconnect();
					}

				} catch (Exception e) {
					e.printStackTrace();
					throw new RuntimeException(e);
				}
			}

		}

	}

}
