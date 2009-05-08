package fi.csc.microarray.messaging.auth;

import javax.swing.SwingUtilities;


public class SimpleAuthenticationRequestListener implements AuthenticationRequestListener {

	private boolean isRequested = false;
	private ClientLoginListener loginListener;
	
	// This class is not used for real usernames and passwords for the users,
	// thus not using char[] etc
	private String username;
	private String password;
	
	public SimpleAuthenticationRequestListener(String username, String password) {
		this.username = username;
		this.password = password;
	}
	
	
	public Credentials authenticationRequest() {
		isRequested = true;
		SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				if (loginListener != null) {
					loginListener.firstLogin();
				}
			}
		});
		return new Credentials(username, password);
	}

	public boolean isRequested() {
		return isRequested;
	}
	
	public void reset() {
		isRequested = false;
	}

	public void setLoginListener(ClientLoginListener loginListener) {
		this.loginListener = loginListener;		
	}

	public void authenticationSucceeded() {
		// ignore	
	}
}
