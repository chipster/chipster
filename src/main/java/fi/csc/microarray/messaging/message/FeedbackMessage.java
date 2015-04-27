package fi.csc.microarray.messaging.message;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import javax.jms.JMSException;
import javax.jms.MapMessage;

import org.apache.log4j.Logger;

import fi.csc.microarray.util.Strings;

/**
 * Message that client sends to manager to give feedback,
 * give details about errors etc.
 * 
 * @author naktinis
 *
 */
public class FeedbackMessage extends ParameterMessage {

    private static final Logger logger = Logger
            .getLogger(JobMessage.class);
    
    private static final String KEY_DETAILS = "details";
    private static final String KEY_EMAIL = "email";
    private static final String KEY_SESSION = "session-url";
    private static final String KEY_LOGS = "logs";
    
    private String details;
    private String email;
    private String url;
    private List<String> logs = new LinkedList<String>();
    
    public FeedbackMessage() {
        super();
    }
    
    public FeedbackMessage(String details, String email, String url) {
        super();
        this.details = details;
        this.email = email;
        this.url = url;
    }
    
    @Override
    public void unmarshal(MapMessage from) throws JMSException {
        super.unmarshal(from);

        // load details
        details = from.getString(KEY_DETAILS);
        email = from.getString(KEY_EMAIL);
        url = from.getString(KEY_SESSION);
        logs = Arrays.asList(!from.getString(KEY_LOGS).equals("") ?
                from.getString(KEY_LOGS).split(";") : new String[] {});
        logger.debug("Unmarshalled " + KEY_DETAILS + " : " + details);
        logger.debug("Unmarshalled " + KEY_EMAIL + " : " + email);
        logger.debug("Unmarshalled " + KEY_SESSION + " : " + url);
        logger.debug("Unmarshalled " + KEY_LOGS + " : " + logs);
    }
    
    @Override
    public void marshal(MapMessage to) throws JMSException {
        super.marshal(to);
        
        // marshal urls into a single string
        String marshalledLog = Strings.delimit(logs, ";");
        
        // log
        logger.debug("Marshalling: " + KEY_DETAILS + " : " + details);
        logger.debug("Marshalling: " + KEY_EMAIL + " : " + email);
        logger.debug("Marshalling: " + KEY_SESSION + " : " + url);
        logger.debug("Marshalling: " + KEY_LOGS + " : " + marshalledLog);
        
        // add details
        to.setString(KEY_DETAILS, details);
        to.setString(KEY_EMAIL, email);
        to.setString(KEY_SESSION, url);
        to.setString(KEY_LOGS, marshalledLog);
    }
    
    public void addLog(String name, String url) {
        logs.add(name + "," + url);
    }
    
    public String getDetails() {
        return details;
    }
    
    public String getEmail() {
        return email;
    }
    
    public String getSessionURL() {
        return url;
    }
    
    /**
     * @return a list of string pairs [(name, url), (name, url), ...].
     */
    public List<String[]> getLogs() {
        List<String[]> logsArrayed = new LinkedList<String[]>();
        for (String log : logs) {
            logsArrayed.add(log.split(","));
        }
        return logsArrayed;
    }
    
}