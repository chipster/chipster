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
    
    private String details;
    private String email;
    
    public FeedbackMessage() {
        super();
    }
    
    public FeedbackMessage(String details, String email) {
        super();
        this.details = details;
        this.email = email;
    }
    
    @Override
    public void unmarshal(MapMessage from) throws JMSException {
        super.unmarshal(from);

        // load details
        this.details = from.getString(KEY_DETAILS);
        this.email = from.getString(KEY_EMAIL);
        logger.debug("Unmarshalled " + KEY_DETAILS + " : " + details);
        logger.debug("Unmarshalled " + KEY_EMAIL + " : " + email);
    }
    
    @Override
    public void marshal(MapMessage to) throws JMSException {
        super.marshal(to);
        
        logger.debug("Marshalling: " + KEY_DETAILS + " : " + details);
        logger.debug("Marshalling: " + KEY_EMAIL + " : " + email);
        
        // add details
        to.setString(KEY_DETAILS, details);
        to.setString(KEY_EMAIL, email);
    }
    
    public String getDetails() {
        return details;
    }
    
    public String getEmail() {
        return email;
    }
    
}