package fi.csc.microarray.filebroker;

/**
 * Java representation of a row in database table 'file' 
 * 
 * @author klemela
 */
public class DbFile {
	
	private String uuid;
	private long size;
	private String created;
	private String lastAccessed;
	
	public DbFile(String uuid, long size, String created, String lastAccessed) {
		this.uuid = uuid;
		this.size = size;
		this.created = created;
		this.lastAccessed = lastAccessed;
	}
	
	public String getUuid() {
		return uuid;
	}
	public void setUuid(String uuid) {
		this.uuid = uuid;
	}
	public long getSize() {
		return size;
	}
	public void setSize(long size) {
		this.size = size;
	}
	public String getCreated() {
		return created;
	}
	public void setCreated(String created) {
		this.created = created;
	}
	public String getLastAccessed() {
		return lastAccessed;
	}
	public void setLastAccessed(String lastAccessed) {
		this.lastAccessed = lastAccessed;
	}
	
	public String toString() {
		return getClass().getSimpleName() + "\tuuid: " + uuid + "\tsize: " + size + "\tcreated: " + created + "\tlastAccessed: " + lastAccessed;
	}
}
