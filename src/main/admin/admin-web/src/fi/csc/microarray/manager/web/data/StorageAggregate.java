package fi.csc.microarray.manager.web.data;

import java.io.Serializable;

import org.springframework.web.util.HtmlUtils;

public class StorageAggregate implements Serializable {
	
	private String username;
	private long size;
	
	public String getUsername() {
		return username;
	}
	public void setUsername(String username) {
		this.username = HtmlUtils.htmlEscape(username);
	}
	public long getSize() {
		return size;
	}
	public void setSize(long size) {
		this.size = size;
	}
}