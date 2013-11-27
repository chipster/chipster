package fi.csc.microarray.filebroker;

/**
 * Java representation of a row in database table 'session'
 * 
 * @author klemela
 */
public class DbSession {
	
	private String uuid;
	private String name;
	private String username;
	
	public DbSession(String uuid, String name, String username) {
		this.uuid = uuid;
		this.name = name;
		this.username = username;
	}
	public String getUuid() {
		return uuid;
	}
	public void setUuid(String uuid) {
		this.uuid = uuid;
	}
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}
	public String getUsername() {
		return username;
	}
	public void setUsername(String username) {
		this.username = username;
	}
	
	/**
	 * Return the name of the server session without any directories. 
	 * 
	 * @param name
	 * @return null if the name is a directory
	 */
	public String getBasename() {
		String basename = name.substring(name.lastIndexOf("/") + 1);
		if (basename.length() > 0) {
			return basename;
		} else {
			return null;
		}
	}
}
