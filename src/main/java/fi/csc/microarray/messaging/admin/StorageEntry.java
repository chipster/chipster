package fi.csc.microarray.messaging.admin;

import java.io.Serializable;
import java.util.Date;

@SuppressWarnings("serial")
public class StorageEntry implements Serializable {
	
	private String id;
	private String username;
	private String name;
	private Date date;
	private long size;
	
	public String getID() {
		return id;
	}
	
	public void setID(String id) {
		this.id = id;
	}
	
	public String getUsername() {
		return username;
	}
	public void setUsername(String username) {
		this.username = username;
	}
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}
	public Date getDate() {
		return date;
	}
	public void setDate(Date date) {
		this.date = date;
	}
	public long getSize() {
		return size;
	}
	public void setSize(long size) {
		this.size = size;
	}
}