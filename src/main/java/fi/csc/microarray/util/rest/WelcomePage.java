package fi.csc.microarray.util.rest;

import java.io.IOException;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServletResponse;

public class WelcomePage {

	private static final String WELCOME_MESSAGE_HEADER = 	
		"<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Strict//EN\"\n" + 
		" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd\">\n" + 
		"<html xmlns=\"http://www.w3.org/1999/xhtml\" xml:lang=\"en\" lang=\"en\">\n" + 
		" <head>\n" + 
		" <meta http-equiv=\"Content-Type\" content=\"text/html; charset=UTF-8\" />\n" + 
		" <title>Chipster</title>\n" + 
		" </head>\n" + 
		" <body>\n" + 
		" <p>Chipster file broker listening at ";
	
	private static final String WELCOME_MESSAGE_FOOTER =
		"</p>" + 
		" </body>\n" + 
		"</html>\n";

	private String url;

	public WelcomePage(String url) {
		this.url = url;
	}

	public void print(HttpServletResponse response) throws ServletException, IOException {
		// write "welcome message"
		response.getWriter().print(WELCOME_MESSAGE_HEADER + url + WELCOME_MESSAGE_FOOTER);
	}
}

