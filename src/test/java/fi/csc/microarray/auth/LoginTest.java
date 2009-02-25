package fi.csc.microarray.auth;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Calendar;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.Random;

import org.testng.Assert;
import org.testng.annotations.Test;

import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.config.ConfigurationLoader.OldConfigurationFormatException;

public class LoginTest {
	
	@Test(groups = {"unit"} )
	public void expirationTest() throws IOException, OldConfigurationFormatException {

		DirectoryLayout.initialiseClientLayout().getConfiguration();			
		
		String validUsername = "valid";
		String expiredUsername = "expired";
		String password = "password";
		
		// create test passwd file and map
		File testPasswdFile = File.createTempFile("chipster-passwd-test", null);
		
		GregorianCalendar tomorrow = new GregorianCalendar();
		tomorrow.set(Calendar.DATE, tomorrow.get(Calendar.DATE) + 1);
		GregorianCalendar yesterday = new GregorianCalendar();
		yesterday.set(Calendar.DATE, yesterday.get(Calendar.DATE) - 1);

		PrintWriter pWriter = new PrintWriter(new BufferedWriter(new FileWriter(testPasswdFile)));
		pWriter.println(validUsername + ":password:" + SimpleFileLoginModule.EXPIRATION_DATE_FORMAT.format(tomorrow.getTime()));
		pWriter.println(expiredUsername + ":password:" + SimpleFileLoginModule.EXPIRATION_DATE_FORMAT.format(yesterday.getTime()));
		pWriter.close();
			
		File testJaasConfigFile = createTestJaasConfig(testPasswdFile);
		AuthenticationProvider authProvider = createAuthProvider(testJaasConfigFile);

		// test
		Assert.assertTrue(authProvider.authenticate(validUsername, password.toCharArray()));
		Assert.assertFalse(authProvider.authenticate(expiredUsername, password.toCharArray()));
	}
	

	/**
	 * Utility method for testing users file anomalies.
	 */
	public void manualTest(File usersFileToTest, String username, String password) throws IOException, OldConfigurationFormatException {

		DirectoryLayout.initialiseClientLayout().getConfiguration();			
		
		File testJaasConfigFile = createTestJaasConfig(usersFileToTest);
		AuthenticationProvider authProvider = createAuthProvider(testJaasConfigFile);

		// test
		Assert.assertTrue(authProvider.authenticate(username, password.toCharArray()));
	}

	private AuthenticationProvider createAuthProvider(File testJaasConfigFile) throws IOException, OldConfigurationFormatException {
		DirectoryLayout.initialiseClientLayout().getConfiguration();			
		AuthenticationProvider authProvider = new JaasAuthenticationProvider();
		System.setProperty("java.security.auth.login.config", testJaasConfigFile.getPath());
		return authProvider;
	}


	private File createTestJaasConfig(File testPasswdFile) throws IOException {
		File testJaasConfigFile = File.createTempFile("chipster-jaas-config-test", null);
		PrintWriter jWriter = new PrintWriter(new BufferedWriter(new FileWriter(testJaasConfigFile)));
		jWriter.println("Chipster {");
	    jWriter.print("fi.csc.microarray.auth.SimpleFileLoginModule sufficient passwdFile=\"");
	    jWriter.print(testPasswdFile.getPath());		
	    jWriter.println("\";");
	    jWriter.println("};");
		jWriter.close();
		return testJaasConfigFile;
	}

	//@Test(groups = {"stress"} )
	// FIXME JAAS seems not to understand that we try to set up different login configs in expirationTest and stressTest
	// (SimpleFileLoginModule gets wrong users file)
	public void stressLoginTest() throws IOException, OldConfigurationFormatException {

		int usernameMaxLength = 5;
		int passwordMaxLength = 5;
		
		int userCount = 200;
		int bruteCount = 100000;
		
		// create test passwd file and map
		HashMap<String, String> users = new HashMap<String, String>();
		File testPasswdFile = File.createTempFile("chipster-passwd-test", null);
		PrintWriter pWriter = new PrintWriter(new BufferedWriter(new FileWriter(testPasswdFile)));
		
		String username;
		String password;
		for (int i = 0; i < userCount; i++) { 
			username = getShortRandomString(usernameMaxLength);
			if (users.get(username) == null) {
				password = getShortRandomString(passwordMaxLength);
				users.put(username, password);
				pWriter.println(username + ':' + password + getDescriptionColumn());
			}
		}
		pWriter.close();
		
		File testJaasConfigFile = createTestJaasConfig(testPasswdFile);
		AuthenticationProvider authProvider = createAuthProvider(testJaasConfigFile);

		// brute force
		int successCount = 0;
		int failCount = 0;
		String bruteUsername;
		String brutePassword;
		String realPassword;
		
		for (int i = 0; i < bruteCount; i++) {
			bruteUsername = getShortRandomString(usernameMaxLength);
			brutePassword = getShortRandomString(passwordMaxLength);

			// auth successful, check that it really is
			if (authProvider.authenticate(bruteUsername, brutePassword.toCharArray())) {
				realPassword = users.get(bruteUsername);
				Assert.assertNotNull(realPassword);
				Assert.assertEquals(realPassword, brutePassword);
				successCount++;
			} 
			
			// auth failed, check that it really should
			else {
				realPassword = users.get(bruteUsername);
				if (realPassword != null) {
					Assert.assertNotSame(realPassword, brutePassword);
				}
				failCount++;
			}
		}
		
		System.out.println("successful: " + successCount + ", failed: " + failCount);
		
	}
	
	
	@Test(groups = {"unit"} )
	public void randomTest() {
		for (int i = 0; i < 100; i++) {
			Assert.assertTrue(getShortRandomString(3).length() <= 3);
		}
	}
	
	private String getShortRandomString(int maxLength) {
		
		Random r = new Random();
		String longString = "";
		String result = "";
		int length;
		boolean valid = false;
		
		while (!valid) {
			longString = randomString();
			length = r.nextInt(10000) % maxLength ;
			result = longString.substring(2, length + 1 + 2);
			if (result.indexOf(':') == -1) {
				valid = true;
			}
		}
		
		return result;
	}
	
	private String randomString() {
		Random r = new Random();
		String token = Long.toString(Math.abs(r.nextLong()), 36);
		return token;
	}
	
	private String getDescriptionColumn() {
		String result = "";
		Random r = new Random();
		
		int choice = r.nextInt(2345) % 3;
		switch (choice) {
		case 0:
			result = "";
			break;
		case 1:
			result = ":";
			break;
		case 2: 
			result = ":" + randomString();
			break;
		}
	
		return result;
	}
}
