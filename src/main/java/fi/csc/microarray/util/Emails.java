package fi.csc.microarray.util;

import java.util.Properties;

import javax.mail.Message;
import javax.mail.MessagingException;
import javax.mail.Session;
import javax.mail.Transport;
import javax.mail.internet.InternetAddress;
import javax.mail.internet.MimeMessage;

public class Emails {
    
    /**
     * Send a regular email.
     * 
     * @param toEmail receiver's email
     * @param replyTo email for reply or null
     * @param subject
     * @param body
     */
    public static void sendEmail(String toEmail, String replyTo, String subject, String body) {
        Session session = Session.getDefaultInstance(new Properties(), null);
        MimeMessage message = new MimeMessage(session);
        try {
            message.addRecipient(Message.RecipientType.TO,
                    new InternetAddress(toEmail));
            message.setSubject(subject);
            message.setText(body);
            
            // set email for reply
            if (replyTo != null) {
                message.addHeader("Reply-To", replyTo);
            }
            
            Transport.send(message);
        } catch (MessagingException e){
            System.err.println("Email could not be sent.");
        }
    }
    
    /**
     * Send a regular email.
     * 
     * @param toEmail receiver's email
     * @param subject
     * @param body
     */
    public static void sendEmail(String toEmail, String subject, String body) {
        sendEmail(toEmail, null, subject, body);
    }
}
