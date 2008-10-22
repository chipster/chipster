/**
 * 
 */
package fi.csc.microarray.messaging.auth;

public interface ClientLoginListener {
	public void firstLogin();
	public void loginCancelled();
}