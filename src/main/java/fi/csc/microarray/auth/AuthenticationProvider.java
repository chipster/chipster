package fi.csc.microarray.auth;

public interface AuthenticationProvider {
	
	public boolean authenticate(String username, char[] password);

}
