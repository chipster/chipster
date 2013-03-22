package fi.csc.microarray.filebroker;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.sql.SQLException;
import java.util.concurrent.TimeUnit;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.apache.commons.io.FileUtils;
import org.eclipse.jetty.servlet.DefaultServlet;
import org.eclipse.jetty.util.IO;
import org.eclipse.jetty.util.URIUtil;
import org.eclipse.jetty.util.log.Log;
import org.eclipse.jetty.util.log.Logger;

import sun.net.www.protocol.http.HttpURLConnection;
import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.util.Files;
import fi.csc.microarray.util.IOUtils;

/**
* <p>Servlet for RESTful file access in Chipster. Extends DefaultServlet and adds support for HTTP PUT and 
* DELETE methods. Also adds Chipster authentication and security checks.</p>
*   
* @author Aleksi Kallio
*
*/
public class RestServlet extends DefaultServlet {

	private Logger logger = Log.getLogger((String)null); // use Jetty default logger
	
	private String cachePath;
	private String publicPath;
	private String storagePath;
	private int cleanUpTriggerLimitPercentage;
	private int cleanUpTargetPercentage;
	private int cleanUpMinimumFileAge;
	
	private String rootUrl;
	private AuthorisedUrlRepository urlRepository;
	private DerbyMetadataServer metadataServer;
	
	

	public RestServlet(String rootUrl, AuthorisedUrlRepository urlRepository, DerbyMetadataServer metadataServer) {
		this.rootUrl = rootUrl;
		this.urlRepository = urlRepository;
		this.metadataServer = metadataServer;
		
		Configuration configuration = DirectoryLayout.getInstance().getConfiguration();
		cachePath = configuration.getString("filebroker", "cache-path");
		storagePath = configuration.getString("filebroker", "storage-path");
		publicPath = configuration.getString("filebroker", "public-path");
		cleanUpTriggerLimitPercentage = configuration.getInt("filebroker", "clean-up-trigger-limit-percentage");
		cleanUpTargetPercentage = configuration.getInt("filebroker", "clean-up-target-percentage");
		cleanUpMinimumFileAge = configuration.getInt("filebroker", "clean-up-minimum-file-age");
	}
	
	@Override
	public void service(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {

		// check request and pass it up to super class
		if (isValidRequest(request)) {
			super.service(request, response);
		} else {
			response.sendError(HttpURLConnection.HTTP_FORBIDDEN);
		}
	};
	
	private boolean isValidRequest(HttpServletRequest request) {

		
		// allow welcome page
		if (isWelcomePage(request)) {
			return true;
		}

		// any other directory is forbidden
		File file = locateFile(request);
		if (file.isDirectory()) {
			return false;
		}

		// allow known request types 
		if (isCacheRequest(request) || isStorageRequest(request) || isPublicRequest(request) ) {
			return true;
		}
		
		// everything else is forbidden
		return false;
	}

	private boolean isWelcomePage(HttpServletRequest request) {
		return "/".equals(request.getPathInfo());
	}

	private boolean isCacheRequest(HttpServletRequest request) {
		return checkRequest(request, cachePath);
	}

	private boolean isStorageRequest(HttpServletRequest request) {
		return checkRequest(request, storagePath);
	}

	private boolean checkRequest(HttpServletRequest request, String prefix) {
		String path = request.getPathInfo();
		
		if (path == null) {
			return false;
		}
		
		if (!path.startsWith("/" + prefix + "/")) {
			return false;
		}

		if (urlRepository.checkFilenameSyntax(path.substring(("/" + prefix + "/").length()))) {
			return true;
		}
		
		return false;
	}

	private boolean isPublicRequest(HttpServletRequest request) {
		String path = request.getPathInfo();
		return (path != null && path.startsWith("/" + publicPath + "/"));
	}

	
	@Override
	protected void doGet(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
		if (logger.isDebugEnabled()) {
			logger.debug("RESTful file access: GET request for " + request.getRequestURI());
		}
		
		// handle welcome page here, delegate rest to super class
		if (isWelcomePage(request)) {
			new WelcomePage(rootUrl).print(response);
			
		} else {			
			
			// touch the file
			File file = locateFile(request);
			file.setLastModified(System.currentTimeMillis());
			
			// touch metadata database
			if (isStorageRequest(request)) {
				String uuid = urlRepository.stripCompressionSuffix(IOUtils.getFilenameWithoutPath(request));
				try {
					metadataServer.markFileAccessed(uuid);
				} catch (SQLException e) {
					throw new ServletException(e);
				}
			}

			// delegate to super class
			super.doGet(request, response);	
		}
		
	}

	
	@Override
	protected void doPut(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
		
		// log
		if (logger.isDebugEnabled()) {
			logger.debug("RESTful file access: PUT request for " + request.getRequestURI());
		}
		
		// check that URL is authorised (authorised URL also implies that quota has been checked)
		if (!urlRepository.isAuthorised(constructUrl(request))) {
			// deny request
			if (logger.isDebugEnabled()) {
				logger.debug("PUT denied for " + constructUrl(request));
			}
			
			response.sendError(HttpURLConnection.HTTP_UNAUTHORIZED);
			return;			
		}
		
		// replace file if it exists
		File file = locateFile(request);
		if (file.exists()) {			
			boolean success = file.delete(); 
			if (!success) {
				response.sendError(HttpURLConnection.HTTP_INTERNAL_ERROR); // file existed and could not be deleted
				return;
			}
		}

		// get file contents
		FileOutputStream out = new FileOutputStream(file);
		InputStream in = request.getInputStream();
		try {
			IO.copy(in, out);
			
		} catch (IOException e) {
			logger.warn(Log.EXCEPTION, e);
			file.delete();
			throw(e);
			
		} finally {
			IOUtils.closeIfPossible(in);
			IOUtils.closeIfPossible(in);
		}
		
		// add file to metadata database, if needed
		if (isStorageRequest(request)) {
			String uuid = urlRepository.stripCompressionSuffix(IOUtils.getFilenameWithoutPath(request));
			long size = file.length();
			try {
				metadataServer.addFile(uuid, size);
				
			} catch (SQLException e) {
				throw new ServletException(e);
			}
		}

		
		// make sure there's enough usable space left after the transfer
		new Thread(new Runnable() {
			@Override
			public void run() {
				try {

					File userDataDir = new File(getServletContext().getRealPath(cachePath));
					long usableSpaceSoftLimit =  (long) ((double)userDataDir.getTotalSpace()*(double)(100-cleanUpTriggerLimitPercentage)/100);

					if (userDataDir.getUsableSpace() <= usableSpaceSoftLimit) {
						logger.info("after put, user data dir usable space soft limit " + FileUtils.byteCountToDisplaySize(usableSpaceSoftLimit) + 
								" (" + (100-cleanUpTriggerLimitPercentage) + "%) reached, cleaning up");
						Files.makeSpaceInDirectoryPercentage(new File(getServletContext().getRealPath(cachePath)), 100-cleanUpTargetPercentage, cleanUpMinimumFileAge, TimeUnit.SECONDS);
						logger.info("after clean up, usable space is: " + FileUtils.byteCountToDisplaySize(new File(getServletContext().getRealPath(cachePath)).getUsableSpace()));
					} 

				} catch (Exception e) {
					logger.warn("could not clean up space after put", e);
				}
			}
		}, "chipster-fileserver-cache-cleanup").start();

		// we return no content
		response.setStatus(HttpURLConnection.HTTP_NO_CONTENT); 
	}
	
	@Override
	protected void doDelete(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
		if (logger.isDebugEnabled()) {
			logger.debug("RESTful file access: DELETE request for " + request.getRequestURI());
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

