package fi.csc.microarray.client.cli;

/**
 * Exception for CliClient for errors that user is most likely able to fix. This interrupts
 * the execution of the current command, but doesn't exit interactive mode.     
 * 
 * @author klemela
 */
public class UserErrorException extends Exception{

	public UserErrorException(String message) {
		super(message);
	}		
}