/*
 * Created on Feb 17, 2005
 *
 */
package fi.csc.microarray.util;

import java.io.IOException;
import java.net.ServerSocket;
import java.net.SocketTimeoutException;

import javax.net.ssl.SSLServerSocket;
import javax.net.ssl.SSLServerSocketFactory;

/**
 * A very usefull class for recreating mysterious SSL initialisation problems.
 * 
 * Run with (replace akallio): -Djavax.net.ssl.keyStore=/home/akallio/.keystore -Djavax.net.ssl.keyStorePassword=microarray
 * 
 * @author akallio
 *
 */
public class SslTest {
	
	public static void main(String[] args) {
		try {
			ServerSocket ss = getSecureSocket(6666);
			ss.setSoTimeout(1); // übershort timeout, we are just testing...
			try {
				ss.accept();
			} catch (SocketTimeoutException e) {
				// this should happen
			}			
			ss.close();
			System.out.println("ok");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	public static SSLServerSocket getSecureSocket(int port) throws IOException {
		SSLServerSocketFactory sslServerSockFactory = (SSLServerSocketFactory)
		SSLServerSocketFactory.
		getDefault();
		return (SSLServerSocket) sslServerSockFactory.createServerSocket(port);
	}
	
}
