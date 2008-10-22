package fi.csc.microarray.messaging.auth;

/**
 * A callback interface for handling authentication requests.
 * @author akallio
 *
 */
public interface AuthenticationRequestListener {
	
	public class Credentials {
		
		public Credentials(String username, String password) {
			this.username = username;
			this.password = password;
		}
		
		public String username;
		public String password;
		
		@Override
		public String toString() {
			return username + " / " + password.replaceAll(".", "*");
		}
	}
	
	public Credentials authenticationRequest();
	
	public void setLoginListener(ClientLoginListener loginListener);

	public void authenticationSucceeded();
}
