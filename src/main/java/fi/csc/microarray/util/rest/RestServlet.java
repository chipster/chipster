package fi.csc.microarray.util.rest;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.mortbay.jetty.servlet.DefaultServlet;
import org.mortbay.log.Log;
import org.mortbay.util.IO;
import org.mortbay.util.URIUtil;

import sun.net.www.protocol.http.HttpURLConnection;

/**
* <p>Adds support for HTTP PUT, MOVE and 
* DELETE methods. If init parameters read-permission-role and write-permission-role
* are defined then all requests are authorized using the defined roles. Also GET methods are
* authorized. </p>
*   
* @author Aleksi Kallio
*
*/
public class RestServlet extends DefaultServlet {

	private File locateFile(HttpServletRequest request) {		
		return new File(getServletContext().getRealPath(URIUtil.addPaths(request.getServletPath(),request.getPathInfo())));		
	}

	@Override
	public void service(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
		File file = locateFile(request);
		if (file.isDirectory()) {
			// directory browsing is not allowed
			response.sendError(HttpURLConnection.HTTP_NOT_FOUND);
		} else {
			super.service(request, response);
		}
	};
	
	@Override
	protected void doGet(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
		if (Log.isDebugEnabled()) {
			Log.debug("RESTful file access: GET request for " + request.getRequestURI());
		}
		
		super.doGet(request, response);
	}
	
	@Override
	protected void doPut(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
		if (Log.isDebugEnabled()) {
			Log.debug("RESTful file access: PUT request for " + request.getRequestURI());
		}
		
		File file = locateFile(request);
		if (file.exists()) {			
			// replace file if it exists
			boolean success = file.delete(); 
			if (!success) {
				response.sendError(HttpURLConnection.HTTP_INTERNAL_ERROR); // file existed and could not be deleted
				return;
			}
		}
		
		FileOutputStream out = new FileOutputStream(file);
		try {
			IO.copy(request.getInputStream(), out);
			
		} catch (IOException e) {
			Log.warn(Log.EXCEPTION, e); // is this obsolete?
			out.close();
			throw(e);
		}
		
		response.setStatus(HttpURLConnection.HTTP_NO_CONTENT); // we return no content
	}
	
	@Override
	protected void doDelete(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
		if (Log.isDebugEnabled()) {
			Log.debug("RESTful file access: DELETE request for " + request.getRequestURI());
		}
		
		File file = locateFile(request);
		
		if (!file.exists()) {
			response.sendError(HttpURLConnection.HTTP_NOT_FOUND); // file not found
			return;
		} 
		
		boolean success = IO.delete(file); // actual delete operation
		
		if (success) {
			response.setStatus(HttpURLConnection.HTTP_NO_CONTENT); // we return no content
		} else {
			response.sendError(HttpURLConnection.HTTP_INTERNAL_ERROR); // could not be deleted due to internal error
		}
	}
	
}

