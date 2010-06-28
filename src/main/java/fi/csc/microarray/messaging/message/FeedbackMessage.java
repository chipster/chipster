package fi.csc.microarray.messaging.message;

import javax.jms.JMSException;
import javax.jms.MapMessage;

import org.apache.log4j.Logger;

/**
 * Message that client sends to manager to give feedback,
 * give details about errors etc.
 * 
 * @author naktinis
 *
 */
public class FeedbackMessage extends PayloadMessage {

    private static final Logger logger = Logger
            .getLogger(JobMessage.class);
    
    private static final String KEY_DETAILS = "details";
    private static final String KEY_EMAIL = "email";
    private static final String KEY_USER = "user";
    
    private String details;
    private String email;
    private String user;
    
    public FeedbackMessage(String details, String email, String user) {
        super();
        this.details = details;
        this.email = email;
        this.user = user;
    }
    
    @Override
    public void unmarshal(MapMessage from) throws JMSException {
        super.unmarshal(from);

        // load details
        this.details = from.getString(KEY_DETAILS);
        this.email = from.getString(KEY_EMAIL);
        this.email = from.getString(KEY_USER);
        logger.debug("Unmarshalled " + KEY_DETAILS + " : " + details);
        logger.debug("Unmarshalled " + KEY_EMAIL + " : " + email);
        logger.debug("Unmarshalled " + KEY_USER + " : " + user);
    }
    
    @Override
    public void marshal(MapMessage to) throws JMSException {
        super.marshal(to);
        
        logger.debug("Marshalling: " + KEY_DETAILS + " : " + details);
        logger.debug("Marshalling: " + KEY_EMAIL + " : " + email);
        logger.debug("Marshalling: " + KEY_USER + " : " + user);
        
        // add details
        to.setString(KEY_DETAILS, details);
        to.setString(KEY_EMAIL, email);
        to.setString(KEY_USER, user);
    }
    
    public String getDetails() {
        return details;
    }
    
    public String getEmail() {
        return email;
    }
    
    public String getUser() {
        return user;
    }
    
}