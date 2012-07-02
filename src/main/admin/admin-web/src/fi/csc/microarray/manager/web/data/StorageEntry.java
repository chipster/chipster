package fi.csc.microarray.manager.web.data;

import java.io.Serializable;

import fi.csc.microarray.messaging.AdminAPI.NodeStatus.Status;

public class StorageEntry implements Serializable {
	
	private String username;
	private String name;
	private Status date;
	private long size;
	
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
	public Status getDate() {
		return date;
	}
	public void setDate(Status date) {
		this.date = date;
	}
	public long getSize() {
		return size;
	}
	public void setSize(long size) {
		this.size = size;
	}
}