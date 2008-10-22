/*
 * Created on Mar 1, 2005
 *
 */
package fi.csc.microarray.client.dialog;

import fi.csc.microarray.client.tasks.Task;

/**
 * @author Aleksi Kallio
 *
 */
public class ErrorDialogUtils {
	
	public static String getMessage(Object param) {
		String msg;
		
		if (param instanceof Exception) {
			msg = ((Exception)param).getLocalizedMessage();
			if (msg == null || "".equals(msg.trim())) {
				// no message => show exception name
				msg = param.getClass().getSimpleName(); 
			}
			
		} else if (param instanceof Task) {
			msg = ((Task)param).getErrorMessage();
			
		} else {
			// if not known type use something else
			msg = null;			
		}
		
		return msg;
	}
	
	public static String getScreenOutput(Object param) {
		if (param instanceof Task) {
			return  ((Task)param).getScreenOutput();
			
		} else {
			return null;			
		}				
	}
}
