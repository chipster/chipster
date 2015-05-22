package fi.csc.microarray.filebroker;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.sql.SQLException;
import java.text.DecimalFormat;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.time.DurationFormatUtils;
import org.eclipse.jetty.servlet.DefaultServlet;
import org.eclipse.jetty.util.IO;
import org.eclipse.jetty.util.URIUtil;
import org.eclipse.jetty.util.log.Log;
import org.eclipse.jetty.util.log.Logger;
import org.joda.time.DateTime;
import org.joda.time.Duration;

import sun.net.www.protocol.http.HttpURLConnection;
import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.filebroker.AuthorisedUrlRepository.Authorisation;
import fi.csc.microarray.util.IOUtils;

/**
* <p>Servlet for RESTful file access in Chipster. Extends DefaultServlet and adds support for HTTP PUT and 
* DELETE methods. Also adds Chipster authentication and security checks.</p>
*   
* @author Aleksi Kallio
*
*/
public class RestServlet extends DefaultServlet {

	private static final String UPLOAD_FILE_EXTENSION = ".upload";

	private Logger logger = Log.getLogger((String)null); // use Jetty default logger
	
	private String cachePath;
	private String publicPath;
	private String storagePath;
	
	private String rootUrl;
	private AuthorisedUrlRepository urlRepository;
	private DerbyMetadataServer metadataServer;
	private DiskCleanUp cacheCleanUp;

	private boolean useChecksums;
	

	// set from configs
	// specify whether get and put requests are logged
	// using jetty debug logging is not very useful as it is so verbose
	private boolean logRest = false;


	public RestServlet(String rootUrl, AuthorisedUrlRepository urlRepository, DerbyMetadataServer metadataServer, DiskCleanUp cacheCleanUp) {
		this.rootUrl = rootUrl;
		this.urlRepository = urlRepository;
		this.metadataServer = metadataServer;
		this.cacheCleanUp = cacheCleanUp;
		
		Configuration configuration = DirectoryLayout.getInstance().getConfiguration();
		cachePath = FileServer.CACHE_PATH;
		storagePath = FileServer.STORAGE_PATH;
		publicPath = configuration.getString("filebroker", "public-path");
		useChecksums = configuration.getBoolean("messaging", "use-checksums");
		if (configuration.getBoolean("filebroker", "log-rest")) {
			logRest = true;
		}
		logger.info("logging rest requests: " + logRest);				
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
	
	private boolean isValidRequest(HttpServletRequest request) throws ServletException {

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

		if (AuthorisedUrlRepository.checkFilenameSyntax(path.substring(("/" + prefix + "/").length()))) {
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
			
			String checksum;
			try {
				checksum = Md5FileUtils.readMd5(file);
				
				if (checksum != null) {
					response.setHeader(ChecksumInputStream.HTTP_CHECKSUM_KEY, checksum);
				}
			} catch (ChecksumParseException e) {
				logger.info("reading checksum file failed", e);
				// continue without checksum
			}
			
			
			// touch metadata database and get checksum
			if (isStorageRequest(request)) {
				String uuid = AuthorisedUrlRepository.stripCompressionSuffix(IOUtils.getFilenameWithoutPath(request));
				try {
					metadataServer.markFileAccessed(uuid);
					
				} catch (SQLException e) {
					throw new ServletException(e);
				}
			}
						
			// delegate to super class
			DateTime before = new DateTime();
			super.doGet(request, response);
			DateTime after = new DateTime();

			// log performance
			Duration duration = new Duration(before, after);
			double rate = getTransferRate(file.length(), duration);
			if (logRest) {
				logger.info("GET " + file.getName()  + " " + 
						"from " + request.getRemoteHost() + " | " +
						FileUtils.byteCountToDisplaySize(file.length()) + " | " + 
						DurationFormatUtils.formatDurationHMS(duration.getMillis()) + " | " +
						new DecimalFormat("###.##").format(rate*8) + " Mbit/s" + " | " +
						new DecimalFormat("###.##").format(rate) + " MB/s");
			}
		}		
	}

	
	@Override
	protected void doPut(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
		
		// log
		if (logger.isDebugEnabled()) {
			logger.debug("RESTful file access: PUT request for " + request.getRequestURI());
		}
		
		// check that URL is authorised (authorised URL also implies that quota has been checked)		
		Authorisation authorisation = urlRepository.getAuthorisation(constructUrl(request));
		if (authorisation == null) {
			// deny request
			if (logger.isDebugEnabled()) {
				logger.debug("PUT denied for " + constructUrl(request));
			}
			
			response.sendError(HttpURLConnection.HTTP_UNAUTHORIZED);
			return;			
		}
		
		File targetFile = locateFile(request);
		
		// create temp file for the upload 
		File tmpFile = getUploadTempFile(targetFile);
		if (tmpFile == null) {
			logger.info("PUT denied for " + constructUrl(request) + ": too many upload temp files for same uuid");
			response.sendError(HttpURLConnection.HTTP_INTERNAL_ERROR);
			return;
		}
		
		//enable one extra byte to recognize misbehaving clients 
		long maxBytes = authorisation.getFileSize() + 1;

		// get file contents
		FileOutputStream out = new FileOutputStream(tmpFile);				
		ChecksumInputStream in = new ChecksumInputStream(request.getInputStream(), useChecksums);				
		
		try {
			DateTime before = new DateTime();
			IO.copy(in, out, maxBytes);
			DateTime after = new DateTime();
			Duration duration = new Duration(before, after);

			double rate = getTransferRate(authorisation.getFileSize(), duration);
			if (logRest) {
				logger.info("PUT " + targetFile.getName()  + " " + 
						"from " + request.getRemoteHost() + " | " +
						FileUtils.byteCountToDisplaySize(authorisation.getFileSize()) + " | " + 
						DurationFormatUtils.formatDurationHMS(duration.getMillis()) + " | " +
						new DecimalFormat("###.##").format(rate*8) + " Mbit/s" + " | " +
						new DecimalFormat("###.##").format(rate) + " MB/s");
			}
			
		} catch (IOException e) {
			logger.warn(Log.EXCEPTION, e);
			tmpFile.delete();
			throw(e);
			
		} finally {
			IOUtils.closeIfPossible(in);
			IOUtils.closeIfPossible(out);
		}
		
		// check that file size matches
		if (tmpFile.length() < authorisation.getFileSize()) {
			logger.info("PUT denied for " + constructUrl(request) + ": stream was shorter than authorised file size");
			response.sendError(HttpURLConnection.HTTP_INTERNAL_ERROR);
			return;
		}
		
		if (tmpFile.length() > authorisation.getFileSize()) {
			logger.info("PUT denied for " + constructUrl(request) + ": stream was longer than authorised file size");
			response.sendError(HttpURLConnection.HTTP_INTERNAL_ERROR);
			return;
		}
		
		String uuid = AuthorisedUrlRepository.stripCompressionSuffix(IOUtils.getFilenameWithoutPath(request));
		long size = tmpFile.length();
		String checksum = in.getChecksum();
		
		if (useChecksums) {
			Md5FileUtils.writeMd5(checksum, targetFile);
		}

		// make file visible		
		if (targetFile.exists()) {
			// someone else was faster to upload the same file, keep it 
			logger.debug("uploaded file exists already, keeping the old one");
			tmpFile.delete();
		} else {
			logger.debug("rename uploaded file to make it visible");
			boolean success = tmpFile.renameTo(targetFile);		
			if (!success) {
				if (targetFile.exists()) {
					// if the rename failed and the file exists, someone else added the file after the exists() call above.
					// in this case the end result is what the client requested, end there is no need to send error
					tmpFile.delete();
				} else {
					logger.debug("rename failed");
					response.sendError(HttpURLConnection.HTTP_INTERNAL_ERROR); // file could not be renamed
					return;
				}
			}
		}
		
		
		// add file to metadata database, if needed
		if (isStorageRequest(request)) {
			try {
				metadataServer.addFile(uuid, size);
			} catch (SQLException e) {
				throw new ServletException(e);
			}
		}
		
		// make sure there's enough usable space left after the transfer
		// (in case of concurrent uploads)
		if (cacheCleanUp != null) {
			cacheCleanUp.spaceRequest(0, false, "Schedule clean up after PUT");
		}

		if (useChecksums) { 
			response.setHeader(ChecksumInputStream.HTTP_CHECKSUM_KEY, checksum);
		}
		// we return no content
		response.setStatus(HttpURLConnection.HTTP_NO_CONTENT); 
	}

	
	private double getTransferRate(long fileSize, Duration duration) {
		double rate;
		if (duration.getMillis() != 0 ) {
			rate = ((double)fileSize) / duration.getMillis()*1000/1024/1024;
		} else {
			rate = 0;
		}

		return rate;
	}


	private File getUploadTempFile(File targetFile) {
		File tmpFile = null;
		for (int i = 1; i < 100; i++) {
			tmpFile = new File(targetFile.getAbsolutePath() + "." + i + UPLOAD_FILE_EXTENSION);
			if (!tmpFile.exists()) {
				break;
			}
		}
		return tmpFile;
	}
	
	@Override
	protected void doDelete(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
		response.setStatus(HttpURLConnection.HTTP_BAD_METHOD);
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
	
	private File locateFile(HttpServletRequest request) throws ServletException {
		String requestPath = URIUtil.addPaths(request.getServletPath(), request.getPathInfo());
		String realPath = getServletContext().getRealPath(requestPath);
		if (realPath == null) {
			throw new ServletException("Servlet context refused to convert a request path to a real path. Make sure that serving files through symlinks is enabled");
		}
		return new File(realPath);		
	}
	
	private URL constructUrl(HttpServletRequest request) throws MalformedURLException {
		return new URL(rootUrl + request.getPathInfo());
	}

}

