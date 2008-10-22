/*
 * Created on Feb 17, 2005
 *
 */
package fi.csc.microarray.util;

import java.io.FileInputStream;
import java.security.Key;
import java.security.KeyStore;
import java.util.Enumeration;

import javax.net.ssl.SSLServerSocketFactory;


/**
 * @author akallio
 *
 */
public class SecurityInfoTool {
	
	private static final String keystoreFile = "/home/akallio/.keystore";
	private static final char[] keystorePass = "microarray".toCharArray();
	
	public static void main(String[] args) {		
		try {
			SSLServerSocketFactory s = (SSLServerSocketFactory)SSLServerSocketFactory.getDefault();
			System.out.print("Default SSL cipher suites: ");
			print(s.getDefaultCipherSuites());
			System.out.print("Supported SSL cipher suites: ");
			print(s.getSupportedCipherSuites());
			
			KeyStore ks = KeyStore.getInstance(KeyStore.getDefaultType());
			ks.load(new FileInputStream(keystoreFile), keystorePass);
			
			System.out.println("Keystore provider info: ");			
			System.out.println(ks.getProvider().getInfo());
			System.out.println("Key aliases in keystore: ");
			Enumeration<String> aliases = ks.aliases();	
			
			while (aliases.hasMoreElements()) {
				String alias = (String)aliases.nextElement();
				Key key = ks.getKey(alias, keystorePass);
				System.out.println(alias + " (algorithm: " + key.getAlgorithm() + ", in provider: " + ks.getProvider().containsKey(key) + ")");
			}
		} catch (Exception e) {
			e.printStackTrace(System.out);
		}
		
	}
	
	private static void print(String [] stuff) {
		for (int i = 0; i < (stuff.length - 1); i++) {
			System.out.print(stuff[i] + ", ");
		}
		System.out.println(stuff[stuff.length-1]);
	}
}
