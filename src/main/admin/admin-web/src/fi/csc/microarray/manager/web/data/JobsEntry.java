package fi.csc.microarray.manager.web.data;

import java.io.Serializable;
import java.util.Date;

import org.springframework.web.util.HtmlUtils;


public class JobsEntry implements Serializable {

	private String operation;
	private String status;
	private Date startTime;
	private String username;
	private String compHost;
		
	public String getOperation() {
		return operation;
	}
	public void setOperation(String operation) {
		this.operation = HtmlUtils.htmlEscape(operation);
	}
	public String getStatus() {
		return status;
	}
	public void setStatus(String status) {
		this.status = HtmlUtils.htmlEscape(status);
	}
	public Date getStartTime() {
		return startTime;
	}
	public void setStartTime(Date startTime) {
		this.startTime = startTime;
	}
	public String getUsername() {
		return username;
	}
	public void setUsername(String username) {
		this.username = HtmlUtils.htmlEscape(username);
	}
	public String getCompHost() {
		return compHost;
	}
	public void setCompHost(String compHost) {
		this.compHost = HtmlUtils.htmlEscape(compHost);
	}
}
