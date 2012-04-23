package fi.csc.microarray.filebroker;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.concurrent.TimeUnit;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.apache.commons.io.IOUtils;
import org.mortbay.jetty.servlet.DefaultServlet;
import org.mortbay.log.Log;
import org.mortbay.util.IO;
import org.mortbay.util.URIUtil;

import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.util.Files;

import sun.net.www.protocol.http.HttpURLConnection;

/**
* <p>Servlet for RESTful file access in Chipster. Extends DefaultServlet and adds support for HTTP PUT and 
* DELETE methods. Also adds Chipster authentication and security checks.</p>
*   
* @author Aleksi Kallio
*
*/
public class RestServlet extends DefaultServlet {

	private String userDataPath;
	private String publicDataPath;
	private int cleanUpFreeSpacePerentage;
	private int cleanUpMinimumFileAge;
	
	private AuthorisedUrlRepository urlRepository;
	private String rootUrl;

	public RestServlet(AuthorisedUrlRepository urlRepository, String rootUrl) {
		this.urlRepository = urlRepository;
		this.rootUrl = rootUrl;
		
		Configuration configuration = DirectoryLayout.getInstance().getConfiguration();
		userDataPath = configuration.getString("filebroker", "user-data-path");
		publicDataPath = configuration.getString("filebroker", "public-data-path");
		cleanUpFreeSpacePerentage = configuration.getInt("filebroker", "clean-up-free-space-percentage");
		cleanUpMinimumFileAge = configuration.getInt("filebroker", "clean-up-minimum-file-age");
	}
	
	@Override
	public void service(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {

		// check request and pass it up to super class
		if (isValidRequest(request)) {
			super.service(request, response);
		} else {
			response.sendError(HttpURLConnection.HTTP_NOT_FOUND);
		}
	};
	
	private boolean isValidRequest(HttpServletRequest request) {

		// must not be directory
		File file = locateFile(request);
		if (file.isDirectory()) {
			return false;
		}

		if (isWelcomePage(request) || isUserDataRequest(request) || isPublicDataRequest(request) ) {
			return true;
		}
		
		return false;
	}

	private boolean isWelcomePage(HttpServletRequest request) {
		return "/".equals(request.getPathInfo());
	}

	private boolean isUserDataRequest(HttpServletRequest request) {
		String path = request.getPathInfo();
		
		if (path == null) {
			return false;
		}
		
		if (!path.startsWith("/" + userDataPath + "/")) {
			return false;
		}

		if (urlRepository.checkFilenameSyntax(path.substring(("/" + userDataPath + "/").length()))) {
			return true;
		}
		
		return false;
	}
	
	private boolean isPublicDataRequest(HttpServletRequest request) {
		String path = request.getPathInfo();
		return (path != null && path.startsWith("/" + publicDataPath + "/"));
	}

	
	@Override
	protected void doGet(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
		if (Log.isDebugEnabled()) {
			Log.debug("RESTful file access: GET request for " + request.getRequestURI());
		}
		
		// handle welcome page here, delegate rest to super class
		if (isWelcomePage(request)) {
			new WelcomePage(rootUrl).print(response);
		} else {
			
			// touch the file
			File file = locateFile(request);
			file.setLastModified(System.currentTimeMillis());
			
			// delegate to super class
			super.doGet(request, response);	
		}
		
	}

	
	@Override
	protected void doPut(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
		if (Log.isDebugEnabled()) {
			Log.debug("RESTful file access: PUT request for " + request.getRequestURI());
		}
		
		// check that URL is authorised
		if (!urlRepository.isAuthorised(constructUrl(request))) {
			// deny request
			if (Log.isDebugEnabled()) {
				Log.debug("PUT denied for " + constructUrl(request));
			}
			
			response.sendError(HttpURLConnection.HTTP_FORBIDDEN);
			return;			
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
		InputStream in = request.getInputStream();
		try {
			IO.copy(in, out);
		} catch (IOException e) {
			Log.warn(Log.EXCEPTION, e); // is this obsolete?
			file.delete();
			throw(e);
		} finally {
			// TODO check that IO.copy and closing
			IOUtils.closeQuietly(out);
			IOUtils.closeQuietly(in);
		}
		
		// make sure there's preferred space left after the transfer
		new Thread(new Runnable() {
			@Override
			public void run() {
				try {
					Files.makeSpaceInDirectoryPercentage(new File(getServletContext().getRealPath(userDataPath)), cleanUpFreeSpacePerentage, cleanUpMinimumFileAge, TimeUnit.SECONDS);
				} catch (Exception e) {
					Log.warn("could not clean up space after put", e);
				}
			}
		}, "chipster-fileserver-cache-cleanup").start();

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
	
	@Override
	protected void doHead(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
		response.setStatus(HttpURLConnection.HTTP_BAD_METHOD);
	};
	
	@Override
	protected void doOptions(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
		response.setStatus(HttpURLConnection.HTTP_BAD_METHOD);
	};
	
	@Override
	protected void doPost(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
		response.setStatus(HttpURLConnection.HTTP_BAD_METHOD);
	};
	
	@Override
	protected void doTrace(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
		response.setStatus(HttpURLConnection.HTTP_BAD_METHOD);
	};
	
	private File locateFile(HttpServletRequest request) {		
		return new File(getServletContext().getRealPath(URIUtil.addPaths(request.getServletPath(), request.getPathInfo())));		
	}
	
	private URL constructUrl(HttpServletRequest request) throws MalformedURLException {
		return new URL(rootUrl + request.getPathInfo());
	}

}

