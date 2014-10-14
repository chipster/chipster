package fi.csc.microarray.messaging.admin;

import java.io.Serializable;
import java.util.Date;

import org.apache.commons.lang3.StringUtils;


public class JobsEntry implements Serializable {

	private String operation;
	private String status;
	private Date startTime;
	private String username;
	private String compHost;
	private String jobId;
		
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
	public void setJobId(String jobId) {
		this.jobId = jobId;
	}
	
	public String getJobId() {
		return jobId;
	}
	
	public String toString() {		
		return ""
				+ StringUtils.rightPad(jobId, 40)
				+ StringUtils.rightPad(compHost, 30)
				+ StringUtils.rightPad(operation, 20)
				+ StringUtils.rightPad(username, 20)
				+ StringUtils.rightPad(startTime.toString(), 30)
				+ StringUtils.rightPad(status, 20);
	}
	public static String getToStringHeaders() {
		return ""
				+ StringUtils.rightPad("JOB ID", 40)
				+ StringUtils.rightPad("HOST", 30)
				+ StringUtils.rightPad("OPERATION", 20)
				+ StringUtils.rightPad("USERNAME", 20)
				+ StringUtils.rightPad("START TIME", 30)
				+ StringUtils.rightPad("STATUS", 20);
	}
}
