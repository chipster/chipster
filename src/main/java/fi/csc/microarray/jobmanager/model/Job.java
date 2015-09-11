package fi.csc.microarray.jobmanager.model;

import java.io.IOException;
import java.util.Date;
import java.util.HashMap;

import javax.jms.JMSException;
import javax.jms.MapMessage;
import javax.persistence.Entity;
import javax.persistence.Id;
import javax.persistence.Lob;

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
public class Job {
	
	public Job() {
		// for Hibernate
	}

	@Id
	private String jobId;
	@Lob
	private String jobMessage;
	@Lob
	private String results;
	private String compId;
	private JobState state;
	private String replyToConnectionId;
	private int replyToSequenceId;

	private Date created;
	private Date scheduled;
	private Date finished;
	private Date seen;

//	private Date rescheduled;
//	private Date dequeued;
//	private Date explicitWait;
//	private String analysisId;
//	private String username;


	Job(JobMessage jobMessage) {
		this.jobId = jobMessage.getJobId();
		setJobMessage(jobMessage);
		setReplyTo((ActiveMQTempTopic) jobMessage.getReplyTo());
		this.created = new Date();
		this.state = JobState.WAITING;
	}

	private void setJobMessage(JobMessage jobMessage) {
		this.jobMessage = toJson(jobMessage);		
	}

	public String getJobId() {
		return jobId;
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

	public ActiveMQTempTopic getReplyTo() {
		ActiveMQTempTopic destination = new ActiveMQTempTopic(new ConnectionId(replyToConnectionId), replyToSequenceId);
		return destination;
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

	void setReplyTo(ActiveMQTempTopic replyTo) {
		this.replyToConnectionId = replyTo.getConnectionId();
		this.replyToSequenceId = replyTo.getSequenceId();
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
			
			if (chipsterMessage.getReplyTo() instanceof ActiveMQTempTopic) {
				ActiveMQTempTopic tempTopic = (ActiveMQTempTopic) chipsterMessage.getReplyTo();
				msgMap.put("replyToConnectionId", tempTopic.getConnectionId());
				msgMap.put("replyToSequenceId", tempTopic.getSequenceId());	
			} else if (chipsterMessage.getReplyTo() instanceof ActiveMQTopic) {
				ActiveMQTopic topic = (ActiveMQTopic) chipsterMessage.getReplyTo();
				msgMap.put("replyToName", topic.getTopicName());	
			} else {
				throw new IllegalArgumentException("unknown type for replyTo " + chipsterMessage.getReplyTo().getClass());
			}
			String json = new GsonBuilder().serializeNulls().create().toJson(msgMap);
			return json;
		} catch (JMSException | IOException e) {
			throw new IllegalArgumentException("unable to marshal chipster message", e);
		}
	}
	
	private ActiveMQMapMessage toMapMessage(String json) {
		ActiveMQMapMessage msg = new ActiveMQMapMessage();
		@SuppressWarnings("unchecked")
		LinkedTreeMap<String, Object> msgMap = new Gson().fromJson(json, LinkedTreeMap.class);
		@SuppressWarnings("unchecked")
		LinkedTreeMap<String, String> properties = (LinkedTreeMap<String, String>) msgMap.get("properties");
		@SuppressWarnings("unchecked")
		LinkedTreeMap<String, String> content = (LinkedTreeMap<String, String>) msgMap.get("content");
		try {
			for (String key : properties.keySet()) {
				msg.setStringProperty(key, properties.get(key));
			}
			for (String key : content.keySet()) {
				msg.setString(key, content.get(key));
			}
			if (msgMap.containsKey("replyToName")) {
				msg.setReplyTo(new ActiveMQTopic((String) msgMap.get("replyToName")));
			} else {
				String connectionId = (String) msgMap.get("replyToConnectionId");
				int sequenceId = (int)(double) msgMap.get("replyToSequenceId");
				msg.setReplyTo(new ActiveMQTempTopic(new ConnectionId(connectionId), sequenceId));
			}
			return msg;
		} catch (JMSException e) {
			throw new IllegalArgumentException("unable to unmarshal chipster message", e);
		}
	}

}
