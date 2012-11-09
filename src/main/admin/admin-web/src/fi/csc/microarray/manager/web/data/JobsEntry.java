package fi.csc.microarray.manager.web.data;

import java.io.Serializable;
import java.util.Date;


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
		this.operation = operation;
	}
	public String getStatus() {
		return status;
	}
	public void setStatus(String status) {
		this.status = status;
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
		this.username = username;
	}
	public String getCompHost() {
		return compHost;
	}
	public void setCompHost(String compHost) {
		this.compHost = compHost;
	}
}
