package fi.csc.microarray.manager.web.data;

import java.io.Serializable;
import java.util.Date;

import javax.persistence.Entity;
import javax.persistence.Id;
import javax.persistence.Table;

@Entity
@Table(name="jobs")
public class JobLogEntry implements Serializable {
	/*
	private static final String CREATE_JOBS_TABLE = 
			"CREATE TABLE IF NOT EXISTS jobs (" +
			"id VARCHAR(100) PRIMARY KEY, " + 
			"operation VARCHAR(200), " +
			"status VARCHAR(200), " + 
			"starttime DATETIME DEFAULT NULL, " + 
			"endtime DATETIME DEFAULT NULL, " +
			"wallclockTime INT DEFAULT NULL, " +
			"errorMessage TEXT DEFAULT NULL, " +
			"outputText TEXT DEFAULT NULL, " + 
			"username VARCHAR(200), " +
			"compHost VARCHAR(500)" +
			");";
		*/
	
	//Id generation method isn't defined, because only database reading is needed
	@Id
	private String id;
	
	//Hibernate takes database column names from field names by default, use @Column(name="EXAMPLE_COLUMN") to override
	private String operation;
	private String state;
	private Date startTime;
	private Date endTime;
	private int wallclockTime;
	private String errorMessage;
	private String outputText;
	private String username;
	private String compHost;
	
	public JobLogEntry() {};
	
	public String getId() {
		return id;
	}
	public void setId(String id) {
		this.id = id;
	}
	public String getOperation() {
		return operation;
	}
	public void setOperation(String operation) {
		this.operation = operation;
	}
	public String getState() {
		return state;
	}
	public void setState(String state) {
		this.state = state;
	}
	public Date getStartTime() {
		return startTime;
	}
	public void setStartTime(Date startTime) {
		this.startTime = startTime;
	}
	public Date getEndTime() {
		return endTime;
	}
	public void setEndTime(Date endTime) {
		this.endTime = endTime;
	}
	public int getWallclockTime() {
		return wallclockTime;
	}
	public void setWallclockTime(int wallclockTime) {
		this.wallclockTime = wallclockTime;
	}
	public String getErrorMessage() {
		return errorMessage;
	}
	public void setErrorMessage(String errorMessage) {
		this.errorMessage = errorMessage;
	}
	public String getOutputText() {
		return outputText;
	}
	public void setOutputText(String outputText) {
		this.outputText = outputText;
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
