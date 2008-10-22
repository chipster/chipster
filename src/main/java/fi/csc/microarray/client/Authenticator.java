package fi.csc.microarray.client;

import java.util.concurrent.CountDownLatch;

import javax.swing.SwingUtilities;

import fi.csc.microarray.client.dialog.LoginDialog;
import fi.csc.microarray.messaging.auth.AuthenticationRequestListener;
import fi.csc.microarray.messaging.auth.ClientLoginListener;

public class Authenticator implements AuthenticationRequestListener {

	private static class CredentialsHolder {
		Credentials credentials = null;
	}
	
	public static interface LoginCallback {
		public boolean authenticate(String username, String password);
	}
		
	private boolean hasBeenAuthenticatedBefore = false;
	private boolean previousAttemptSuccessful = true;
	
	private ClientLoginListener loginListener = null;
	
	
	public void setLoginListener(ClientLoginListener listener) {
		this.loginListener = listener;
	}

	public Credentials authenticationRequest() {
		final CredentialsHolder credentialsHolder = new CredentialsHolder();
		final CountDownLatch loginLatch = new CountDownLatch(1);
		LoginCallback callback = new LoginCallback() {

			public boolean authenticate(String username, String password) {
				// store credentials 
				credentialsHolder.credentials = new Credentials(username, new String(password));
				
				// notify blocking caller
				loginLatch.countDown();
							
				return true;
			}
		};

		if (previousAttemptSuccessful) {
			previousAttemptSuccessful = false;
			new LoginDialog(callback).setVisible(true);
		} else {
			new LoginDialog(callback, true).setVisible(true);
		}
		
		try {
			loginLatch.await();
		} catch (InterruptedException e) {			
			throw new RuntimeException(e);
		}
		
		return credentialsHolder.credentials;
	}

	public void authenticationSucceeded() {		
		// notify blocking UI initialisation, if exists
		if (!hasBeenAuthenticatedBefore) {
			SwingUtilities.invokeLater(new Runnable() {
				public void run() {
					loginListener.firstLogin();
				}
			});
			hasBeenAuthenticatedBefore = true;
		}
		previousAttemptSuccessful = true;
		
	}

}
