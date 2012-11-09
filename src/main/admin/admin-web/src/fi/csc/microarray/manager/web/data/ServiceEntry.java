package fi.csc.microarray.manager.web.data;

import java.io.Serializable;

import org.springframework.web.util.HtmlUtils;

import fi.csc.microarray.messaging.AdminAPI.NodeStatus.Status;

public class ServiceEntry implements Serializable {
	
	private String name;
	private String host;
	private Status status;
	private int count = 0;
	
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = HtmlUtils.htmlEscape(name);
	}
	public String getHost() {
		return host;
	}
	public void setHost(String host) {
		this.host = HtmlUtils.htmlEscape(host);
	}
	public Status getStatus() {
		return status;
	}
	public void setStatus(Status status) {
		this.status = status;
	}
	public int getCount() {
		return count;
	}
	public void setCount(int count) {
		this.count = count;
	} 
}