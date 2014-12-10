package fi.csc.microarray.messaging.message;

import java.util.HashMap;

import javax.jms.JMSException;
import javax.jms.MapMessage;

import com.google.gson.Gson;
import com.google.gson.JsonObject;

public class JsonMessage extends ChipsterMessage {
	
	private final static String KEY_JSON = "json";
	public static final String KEY_STATUS_REPORT = "status-report";
	
	private String value;
	
	public JsonMessage() {
		super();
	}
	
	public JsonMessage(String value) {
		super();
		this.value = value;
	}
	
	public String getJson() {
		return value;
	}
	
	public void unmarshal(MapMessage from) throws JMSException {
		super.unmarshal(from);
		this.value = from.getString(KEY_JSON);
	}

	public void marshal(MapMessage mapMessage) throws JMSException {
		super.marshal(mapMessage);
		mapMessage.setString(KEY_JSON, this.value);
	}

	public void setCommnad(String command) {
		HashMap<String, String> jsonMap = new HashMap<>();
		jsonMap.put("command", command);		
		this.value = new Gson().toJson(jsonMap);
	}
	
	public String getValue(String key) {
		JsonObject jobj = new Gson().fromJson(this.value, JsonObject.class);
		return jobj.get(key).toString();
	}
}