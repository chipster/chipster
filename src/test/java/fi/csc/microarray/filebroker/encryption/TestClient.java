package fi.csc.microarray.filebroker.encryption;

import java.io.ByteArrayInputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.HttpURLConnection;
import java.net.URL;
import java.security.InvalidKeyException;
import java.security.NoSuchAlgorithmException;
import java.util.zip.Deflater;
import java.util.zip.DeflaterOutputStream;

import javax.crypto.Cipher;
import javax.crypto.CipherOutputStream;
import javax.crypto.KeyGenerator;
import javax.crypto.NoSuchPaddingException;
import javax.crypto.SecretKey;
import javax.crypto.spec.SecretKeySpec;
import javax.net.ssl.HostnameVerifier;
import javax.net.ssl.HttpsURLConnection;
import javax.net.ssl.SSLSession;

public class TestClient {

    public static void main(String... args) throws Exception {
        final File file = new File("test_input");
        DataInputStream dis = new DataInputStream(new FileInputStream(file));
        byte[] bytes = new byte[(int) file.length()];
        dis.readFully(bytes);
        
        System.out.println("Test started...");
        burnAndTest(bytes, null, false);
//        burnAndTest(bytes, null, true);
        burnAndTest(bytes, "SSL", false);
//      burnAndTest(bytes, "SSL", true);
        burnAndTest(bytes, "RC4", false);
//        burnAndTest(bytes, "RC4", true);
        burnAndTest(bytes, "AES", false);
//        burnAndTest(bytes, "AES", true);
    }

    private static void burnAndTest(byte[] bytes, String encryptionAlgorithm, boolean compression) throws IOException, InvalidKeyException, NoSuchAlgorithmException, NoSuchPaddingException {
    	
    	// burn in
    	test(bytes, encryptionAlgorithm, compression);
    	
    	// test
    	int runCount = 5;
    	float bandwidthSum = 0f;
    	for (int i = 0; i < runCount; i++) {
    		bandwidthSum += test(bytes, encryptionAlgorithm, compression);
    	}
    	float bandwidth = Math.round(bandwidthSum / (float)runCount);
    	
		System.out.println(encryptionAlgorithm + "/" + (compression ? "compressed" : "normal") + ":  Bandwidth " + bandwidth + " MB/s.");

    }
    
    private static float test(byte[] bytes, String encryptionAlgorithm, boolean compression) throws IOException, InvalidKeyException, NoSuchAlgorithmException, NoSuchPaddingException {

    	String host = "localhost";
    	
    	HttpURLConnection connection;
    	if ("SSL".equals(encryptionAlgorithm)) {
        	System.setProperty("javax.net.ssl.trustStore", "keystore.ks");
        	System.setProperty("javax.net.ssl.trustStorePassword", "microarray");
        	connection = (HttpsURLConnection)new URL("https://" + host + ":9443/test_output_" + encryptionAlgorithm + "_" + compression).openConnection();
        	((HttpsURLConnection)connection).setHostnameVerifier(new HostnameVerifier() {
    			@Override
    			public boolean verify(String hostname, SSLSession session) {
    				return true;
    			}
        	});
        	
    	} else {
        	connection = (HttpURLConnection)new URL("http://" + host + ":9080/test_output_" + encryptionAlgorithm + "_" + compression).openConnection();
    	}
    	
    	connection.setDoOutput(true);
    	connection.setRequestMethod("PUT");
    	connection.setChunkedStreamingMode(2048);
    	
        OutputStream out = connection.getOutputStream();
    	InputStream in = new ByteArrayInputStream(bytes);
        
        long start = System.nanoTime();

        //
        // ENCRYPT
        //
        
        // Get the KeyGenerator
        if (encryptionAlgorithm != null && !"SSL".equals(encryptionAlgorithm)) {
        	KeyGenerator kgen = KeyGenerator.getInstance(encryptionAlgorithm);
        	kgen.init(128); // 192 and 256 bits may not be available

        	// Generate the secret key specs.
        	SecretKey skey = kgen.generateKey();
        	byte[] raw = skey.getEncoded();
        	SecretKeySpec skeySpec = new SecretKeySpec(raw, encryptionAlgorithm);

        	// Instantiate the cipher
        	Cipher cipher = Cipher.getInstance(encryptionAlgorithm);
        	cipher.init(Cipher.ENCRYPT_MODE, skeySpec);

        	// Bytes written to out will be encrypted
        	out = new CipherOutputStream(out, cipher);
        }


        //
        // COMPRESS
        //
        if (compression) {
        	Deflater def = new Deflater(Deflater.BEST_SPEED);
        	out = new DeflaterOutputStream(out, def, 4 * 1024);
        }

        //
        // COPY
        //
        
		// Read in the cleartext bytes and write to out to encrypt
		int numRead = 0;
		int byteCount = 0;
		byte[] buf = new byte[1024];
		while ((numRead = in.read(buf)) >= 0) {
			out.write(buf, 0, numRead);
			byteCount += numRead;
		}
		out.flush();
		out.close();
        
		//
		// OUTPUT RESULTS
		//
		float time = (float)(System.nanoTime() - start);
		float bandwidth = (byteCount * 1000f) / time; // 1 byte per ns = 1000 MB/s.
		return bandwidth; 
    }
    
}
