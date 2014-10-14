package fi.csc.microarray.messaging.admin;

import java.io.Serializable;

public class StorageAggregate implements Serializable {
	
	private String username;
	private long size;
	
	public String getUsername() {
		return username;
	}
	public void setUsername(String username) {
		this.username = username;
	}
	public long getSize() {
		return size;
	}
	public void setSize(long size) {
		this.size = size;
	}
}