package fi.csc.microarray.filebroker;

/**
 * Java representation of a row in database table 'session'
 * 
 * @author klemela
 */
public class DbSession {
	
	private String dataId;
	private String name;
	private String username;
	
	public DbSession(String sessionId, String name, String username) {
		this.dataId = sessionId;
		this.name = name;
		this.username = username;
	}
	public String getDataId() {
		return dataId;
	}
	public void setDataId(String uuid) {
		this.dataId = uuid;
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
