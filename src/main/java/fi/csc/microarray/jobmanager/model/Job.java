package fi.csc.microarray.jobmanager.model;

import java.io.IOException;
import java.util.Date;
import java.util.HashMap;
import java.util.UUID;

import javax.jms.Destination;
import javax.jms.JMSException;
import javax.jms.MapMessage;
import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.EnumType;
import javax.persistence.Enumerated;
import javax.persistence.Id;
import javax.persistence.Index;
import javax.persistence.Lob;
import javax.persistence.Table;

import org.apache.activemq.command.ActiveMQMapMessage;
import org.apache.activemq.command.ActiveMQTempTopic;
import org.apache.activemq.command.ActiveMQTopic;
import org.apache.activemq.command.ConnectionId;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.internal.LinkedTreeMap;

import fi.csc.microarray.messaging.JobState;
import fi.csc.microarray.messaging.message.ChipsterMessage;
import fi.csc.microarray.messaging.message.JobMessage;
import fi.csc.microarray.messaging.message.ResultMessage;

@Entity
@Table(indexes = @Index(name="created", columnList = "created"))
public class Job {
	
	public Job() {
		// for Hibernate
	}

	@Id
	@Column( columnDefinition = "uuid", updatable = false ) // uuid instead of binary
	private UUID jobId;
	@Lob
	private String jobMessage;
	@Lob
	private String results;
	private String compId;
	@Enumerated(EnumType.STRING) // string instead of int
	private JobState state;

	private Date created;
	private Date scheduled;
	private Date finished;
	@SuppressWarnings("unused")
	private Date seen;
	private String replyToConnectionId;
	private int replyToSequenceId;
	private String replyToName;
	private String compHost;

//	private Date rescheduled;
//	private Date dequeued;
//	private Date explicitWait;
//	private String analysisId;
//	private String username;


	Job(JobMessage jobMessage) {
		this.jobId = UUID.fromString(jobMessage.getJobId());
		setJobMessage(jobMessage);
		setReplyTo((ActiveMQTempTopic) jobMessage.getReplyTo());
		this.created = new Date();
		this.state = JobState.WAITING;
	}

	private void setJobMessage(JobMessage jobMessage) {
		this.jobMessage = toJson(jobMessage);
	}

	public String getJobId() {
		return jobId.toString();
	}

	public JobMessage getJobMessage() {
		MapMessage mapMessage = toMapMessage(jobMessage);
		JobMessage jobMessage = new JobMessage();
		try {
			jobMessage.unmarshal(mapMessage);
			return jobMessage;
		} catch (JMSException e) {
			throw new IllegalArgumentException("unable to unmarshal chipster message", e);
		}
	}

	public Date getScheduled() {
		return scheduled;
	}

	public ResultMessage getResults() {				
		MapMessage mapMessage = toMapMessage(results);
		ResultMessage resultMessage = new ResultMessage();
		try {
			resultMessage.unmarshal(mapMessage);
			return resultMessage;
		} catch (JMSException e) {
			throw new IllegalArgumentException("unable to unmarshal chipster message", e);
		}
	}

	public JobState getState() {
		return state;
	}

	public String getCompId() {
		return compId;
	}

	public Destination getReplyTo() {
		if (this.replyToName != null) {
			return new ActiveMQTopic(this.replyToName);
		} else {
			return new ActiveMQTempTopic(new ConnectionId(this.replyToConnectionId), replyToSequenceId);
		}
	}

	public Date getCreated() {
		return created;
	}

	public Date getFinished() {
		return finished;
	}

	public void setFinished(Date finished) {
		this.finished = finished;
	}

	public long getSecondsSinceCreated() {
		return (System.currentTimeMillis() - created.getTime()) / 1000;
	}

	void setScheduled(Date scheduled) {
		this.scheduled = scheduled;
	}

	void setResults(ResultMessage results) {		
		this.results = toJson(results);
	}

	void setState(JobState state) {
		this.state = state;
	}

	void setCompId(String compId) {
		this.compId = compId;
	}

	void setReplyTo(Destination replyTo) {
		if (replyTo instanceof ActiveMQTempTopic) {
			ActiveMQTempTopic tempTopic = (ActiveMQTempTopic) replyTo;
			this.replyToConnectionId = tempTopic.getConnectionId();
			this.replyToSequenceId = tempTopic.getSequenceId();
			this.replyToName = null;
		} else if (replyTo instanceof ActiveMQTopic) {
			ActiveMQTopic topic = (ActiveMQTopic) replyTo;
			try {
				this.replyToConnectionId = null;
				this.replyToSequenceId = 0;
				this.replyToName = topic.getTopicName();
			} catch (JMSException e) {
				throw new IllegalArgumentException("unable to get topic name", e);
			}
		} else {
			throw new IllegalArgumentException("unknown destination type " + replyTo.getClass().getName());
		}
	}

	void setSeen(Date seen) {
		this.seen = seen;
	}

	private String toJson(ChipsterMessage chipsterMessage) {
		ActiveMQMapMessage msg = new ActiveMQMapMessage();
		try {
			chipsterMessage.marshal(msg);
			HashMap<String, Object> msgMap = new HashMap<>();
			HashMap<String, String> properties = new HashMap<>();
			HashMap<String, String> content = new HashMap<>();
			for (String name : msg.getProperties().keySet()) {
				properties.put(name, msg.getStringProperty(name));
			}
			for (String name : msg.getContentMap().keySet()) {
				content.put(name, msg.getString(name));
			}
			
			msgMap.put("properties", properties);
			msgMap.put("content", content);
						
			String json = new GsonBuilder().serializeNulls().create().toJson(msgMap);
			return json;
		} catch (JMSException | IOException e) {
			throw new IllegalArgumentException("unable to marshal chipster message", e);
		}
	}
	
	@SuppressWarnings("unchecked")
	private ActiveMQMapMessage toMapMessage(String json) {
		ActiveMQMapMessage msg = new ActiveMQMapMessage();
		LinkedTreeMap<String, Object> msgMap = new Gson().fromJson(json, LinkedTreeMap.class);
		LinkedTreeMap<String, String> properties = (LinkedTreeMap<String, String>) msgMap.get("properties");
		LinkedTreeMap<String, String> content = (LinkedTreeMap<String, String>) msgMap.get("content");
		try {
			for (String key : properties.keySet()) {
				msg.setStringProperty(key, properties.get(key));
			}
			for (String key : content.keySet()) {
				msg.setString(key, content.get(key));
			}
			return msg;
		} catch (JMSException e) {
			throw new IllegalArgumentException("unable to unmarshal chipster message", e);
		}
	}

	public void setCompHost(String compHost) {
		this.compHost = compHost;
	}
	
	public String getCompHost() {
		return this.compHost;
	}
}
