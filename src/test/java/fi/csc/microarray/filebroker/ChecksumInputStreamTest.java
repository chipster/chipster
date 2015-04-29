package fi.csc.microarray.filebroker;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.URLConnection;
import java.util.Random;

import javax.xml.bind.DatatypeConverter;

import org.junit.Assert;
import org.junit.Test;

import fi.csc.microarray.util.IOUtils;

public class ChecksumInputStreamTest {
	
	String[][] examples = new String[][] {
			{"9E204EB216C2ECBAE690BFFB004E81B2CAA7DB42D6E5773D6C026ABAB65025900614B5A65C56BC0BE08D0D21ABA34D7381A0", "32c027190a78f0fb2ae6819b77cc2864" },
			//check that leading zeros in md5 are preserved
			{"31DB475415008A0D5DF7335E81594F7015E2D546AEEFAC67308C04C068D01325371712A717D83252151780C6692688BA8B4A", "00fdcd1c90cf409169750ea88ad22899" },
			{"C15E04809BF8DBC09395EB7C42CF80B0C99F8B50CFD4426D35FE32398CEF94B009F7C52DC7F4BF9B51EAE7ECB05957115BF1", "04cc1d49414deeb72cbe182cc412de25" }
	};
	
	@Test
	public void test() throws IOException, ChecksumException, ContentLengthException {
		
		for (String[] example : examples) {
			// example md5 is calculated from bytes, not from hex String, 
			// but the hex format is needed to present those bytes in this source code
			byte[] data = DatatypeConverter.parseHexBinary(example[0]);
			
			// test md5 calculation
			String md5 = testCalculation(data);			
			Assert.assertEquals(example[1], md5);
			
			// successful verification
			testVerification(data, example[1]);
			Assert.assertEquals(example[1], md5);
			
			// failing verification
			try {
				testVerification(data, "wrong-md5");
				Assert.fail();
			} catch (ChecksumException e) {				
			}
		}
	}
	
	private String testCalculation(byte[] data) throws IOException {
		
		ChecksumInputStream stream = preprocess(data, null);		
		String md5 =  stream.getChecksum();		
		postprocess(stream);
		
		return md5;
	}
	
	private String testVerification(byte[] data, String exampleMd5) throws IOException, ChecksumException, ContentLengthException {
		
		ChecksumInputStream stream = preprocess(data, new TestURLConnection(exampleMd5));		
		String md5 =  stream.verifyChecksums();		
		postprocess(stream);
		
		return md5;
	}

	private void postprocess(ChecksumInputStream stream) throws IOException {
		try {
			stream.read();
			Assert.fail();
		} catch (IllegalStateException e) {			
		}
		
		try {
			stream.read(new byte[10], 0, 10);
			Assert.fail();
		} catch (IllegalStateException e) {			
		}
		
		try {
			stream.getChecksum();
			Assert.fail();
		} catch (IllegalStateException e) {			
		}
	}

	private ChecksumInputStream preprocess(byte[] data, URLConnection connection) throws IOException {
		ByteArrayInputStream baseStream = new ByteArrayInputStream(data);
		
		ChecksumInputStream stream = new ChecksumInputStream(baseStream, true, connection);
		
		Assert.assertEquals(baseStream.available(), stream.available());
		
		byte[] streamBytes = new byte[data.length];		
		for (int i = 0; i < data.length; i++) {
			streamBytes[i] = (byte) stream.read();
		}
		
		Assert.assertEquals(DatatypeConverter.printHexBinary(data), DatatypeConverter.printHexBinary(streamBytes));		
		Assert.assertEquals(baseStream.available(), stream.available());
		
		stream.close();
		return stream;
	}
	
	public static class TestURLConnection extends URLConnection {

		private String md5;

		public TestURLConnection(String md5) {
			super(null);
			this.md5 = md5;
		}

		@Override
		public void connect() throws IOException {
		}
		
		@Override
		public String getHeaderField(String name) {
	        Assert.assertEquals("Etag", name);
	        return md5;
	    }
	}
	
	/**
	 * Make a quick performance test
	 * 
	 * Example result on ThinkPad X220 with 512 MB of data,
	 * 'Checksum enabled' was constrained by CPU single core performance):
	 * 
	 * Base stream: 	47 ms, 	2127.6597 MB/s
	 * Checksum disabled: 	50 ms, 	2000.0 MB/s
	 * Checksum enabled: 	486 ms, 	205.76132 MB/s
	 * 
	 * @param args
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public static void main(String args[]) throws IOException, InterruptedException {
		
		Random rand = new Random();
		//requires JVM argument -Xmx2g
//		byte[] data = new byte[1024*1024*512];		
		byte[] data = new byte[1024*1024*100];
		rand.nextBytes(data);
		
		ByteArrayInputStream baseStream = null;
		ChecksumInputStream inStream = null;
		ByteArrayOutputStream outStream = null;
		
		baseStream = new ByteArrayInputStream(data);
		outStream = new ByteArrayOutputStream(data.length);		
		measure("Base stream", data.length, baseStream, outStream);
		
		baseStream = null;
		inStream = null;
		outStream = null;
		System.gc();
		
		// sleep for a while to see bottlenecks in system monitor (use big data)
		Thread.sleep(1000);
		
		baseStream = new ByteArrayInputStream(data);
		inStream = new ChecksumInputStream(baseStream, false);
		outStream = new ByteArrayOutputStream(data.length);		
		measure("Checksum disabled", data.length, inStream, outStream);
		
		baseStream = null;
		inStream = null;
		outStream = null;
		System.gc();
		
		Thread.sleep(1000);
		
		baseStream = new ByteArrayInputStream(data);
		inStream = new ChecksumInputStream(baseStream, true);
		outStream = new ByteArrayOutputStream(data.length);		
		measure("Checksum enabled", data.length, inStream, outStream);		
	}

	private static void measure(String message, long bytes, InputStream in,
			OutputStream out) throws IOException {
		long t = System.currentTimeMillis();		
		IOUtils.copy(in, out);
		long dt = System.currentTimeMillis() - t;
		in.close();
		out.close();
		System.out.println(message + ": \t" + dt + " ms, \t" + bytes / 1024 / 1024 / (dt / 1000f) + " MB/s");
	}
}
