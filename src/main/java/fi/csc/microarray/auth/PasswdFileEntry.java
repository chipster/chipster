package fi.csc.microarray.auth;

public class PasswdFileEntry {

	private String username;
	private char[] password;
	
	public PasswdFileEntry(String username, char[] password) {
		this.username = username;
		this.password = password;
	}

	public String getUsername() {
		return username;
	}

	public char[] getPassword() {
		return password;
	}
	
	public void clear() {
		username = null;
		if (password != null) {
			for (int i = 0; i < password.length; i++) {
				password[i] = ' ';
			}
		}
	}
}
