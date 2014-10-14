package fi.csc.chipster.web.adminweb.util;

import com.vaadin.server.Page;
import com.vaadin.ui.Notification;
import com.vaadin.ui.Notification.Type;

import fi.csc.microarray.messaging.message.SuccessMessage;

public class Notificationutil {
	public static void showFailNotification(String title, String description) {
		Notification notification = new Notification(title + "\n", description, Type.WARNING_MESSAGE);
		notification.setDelayMsec(-1);
		notification.setHtmlContentAllowed(false);
		notification.show(Page.getCurrent());
	}
	
	public static void showFailNotification(String title, SuccessMessage message) {
		String description = "";
		String lineBreak = "\n\n";
		if (message.getErrorMessage() != null && !message.getErrorMessage().isEmpty()) {
			description += message.getErrorMessage() + lineBreak;
		}
		
		if (message.getDetails() != null && !message.getDetails().isEmpty()) {
			description += message.getDetails() + lineBreak;
		}

		if (message.getExceptionString() != null && !message.getExceptionString().isEmpty()) {
			description += message.getExceptionString() + lineBreak;
		}
		
		if (description.endsWith(lineBreak)) {
			description = description.substring(0, description.length() - lineBreak.length());
		}
		showFailNotification(title, description);
	}

}
