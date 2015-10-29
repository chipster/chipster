package fi.csc.chipster.tools.common;

import java.io.File;
import java.io.InputStream;
import java.net.URL;
import java.net.URLConnection;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;

import org.apache.commons.io.FileUtils;

import fi.csc.microarray.comp.java.JavaCompJobBase;
import fi.csc.microarray.messaging.JobState;
import fi.csc.microarray.util.Exceptions;
import fi.csc.microarray.util.ToolUtils;
import fi.csc.microarray.util.KeyAndTrustManager;
import fi.csc.microarray.util.UrlTransferUtil;

public class DownloadFile extends JavaCompJobBase {
	
	public static final int CONNECTION_TIMEOUT = 15*1000; //ms
	public static final int READ_TIMEOUT = 7*24*60*60*1000; //ms
	
	/*
	 * We definitely don't wan't to enable 'file' protocol. Other
	 * protocols might be fine as long as the isLocalhost() check below
	 * is sufficient for it.
	 */			
	public static final List<String> allowedProtocols = Arrays.asList(new String[] {"http", "https", "ftp"});
	
	public static final String CURRENT = "current";
	public static final String YES = "yes";
	public static final String NO = "no";
	
	@Override
	public String getSADL() {
		return 	"TOOL DownloadFile.java: \"Download file from URL directly to server\" (Download a file from a URL address to the Chipster server. The URL must be visible to Chipster server. If it's not, use client's 'Import from URL' functionality instead.)" + "\n" +
				"OUTPUT downloaded_file: \"Downloaded file\"\n" +
				"PARAMETER paramUrl: \"URL\" TYPE STRING (URL to download)\n" + 
				"PARAMETER OPTIONAL paramFileExtension: \"Add a file extension\" TYPE [" + CURRENT + ": \"Keep current\", bam: \"BAM\", fa: \"FASTA\", fastq: \"FASTQ\", gtf: \"GTF\"] DEFAULT " + CURRENT + " (The output file is named according to the last part of the URL. If it doesn't contain a correct file extension, select it here so that the file type is recognized correctly.)\n" + 
				"PARAMETER OPTIONAL paramCheckCerts: \"Require valid SSL certificate\" TYPE [" + YES + ": \"Yes\", " + NO + ": \"No\"] DEFAULT " + YES + " (Disable if the https server has a self-signed ssl certificate.)\n";	
	}
	
	@Override
	protected void execute() { 
		updateStateToClient(JobState.RUNNING, "downloading");

		try {
			// file 
			File outputFile = new File(jobWorkDir, toolDescription.getOutputFiles().get(0).getFileName().getID()); 

			// parameters
			List<String> parameters = inputMessage.getParameters(JAVA_PARAMETER_SECURITY_POLICY, toolDescription);
			String urlString = parameters.get(0);
			String fileExtension = parameters.get(1);
			boolean checkCerts = true;
			if (NO.equals(parameters.get(2))) {
				checkCerts = false;
			}
			
			URL url = new URL(urlString);
			
			
			if (!allowedProtocols.contains(url.getProtocol())) {
				
				getResultMessage().setErrorMessage("Unsupported protocol: " + url.getProtocol());
				updateState(JobState.FAILED_USER_ERROR, "");
				return;
			}
			
			if (UrlTransferUtil.isLocalhost(url.getHost())) {

				getResultMessage().setErrorMessage("Not allowed to connect localhost: " + url.getHost());
				updateState(JobState.FAILED_USER_ERROR, "");
				return;
			}
			
			String datasetName = UrlTransferUtil.parseFilename(url);
			
			if (!CURRENT.equals(fileExtension)) {
				datasetName += "." + fileExtension;
			}
			
			URLConnection connection = url.openConnection();
			if (checkCerts) {
				KeyAndTrustManager.configureForCACertificates(connection);
			} else {
				KeyAndTrustManager.configureForTrustAllCertificates(connection);
			}
			connection.setConnectTimeout(CONNECTION_TIMEOUT);
	        connection.setReadTimeout(READ_TIMEOUT);
			InputStream stream = connection.getInputStream();
			FileUtils.copyInputStreamToFile(stream, outputFile);
			
			LinkedHashMap<String, String> nameMap = new LinkedHashMap<>();
			nameMap.put(outputFile.getName(), datasetName);
			
			ToolUtils.writeOutputDescription(jobWorkDir, nameMap);
						
		} catch (Exception e) {
			getResultMessage().setErrorMessage(Exceptions.getStackTrace(e));
			updateState(JobState.FAILED, "");
			return;
		}

		updateStateToClient(JobState.RUNNING, "download finished");
	}
}
