package fi.csc.chipster.tools.common;

import java.io.File;
import java.net.URL;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;

import org.apache.commons.io.FileUtils;

import fi.csc.microarray.analyser.java.JavaAnalysisJobBase;
import fi.csc.microarray.messaging.JobState;
import fi.csc.microarray.util.Exceptions;
import fi.csc.microarray.util.JavaToolUtils;
import fi.csc.microarray.util.UrlTransferUtil;

public class DownloadFile extends JavaAnalysisJobBase {
	
	public static final int CONNECTION_TIMEOUT = 60*1000; //ms
	public static final int READ_TIMEOUT = 7*24*60*60*1000; //ms
	
	/*
	 * We definitely don't wan't to enable 'file' protocol. Other
	 * protocols might be fine as long as the isLocalhost() check below
	 * is sufficient for it.
	 */			
	public static final List<String> allowedProtocols = Arrays.asList(new String[] {"http", "https", "ftp"});
	
	public static final String CURRENT = "current";
	
	@Override
	public String getSADL() {
		return 	"TOOL DownloadFile.java: \"Download file\" (Download a file from an URL address to the Chipster server. The URL must be visible to Chipster server. If it's not, use client's 'Import from URL' functionality instead.)" + "\n" +
				"OUTPUT downloaded_file: \"Downloaded file\"" + "\n" +
				"PARAMETER paramUrl: \"URL\" TYPE STRING (URL to download)" + 
				"PARAMETER paramFileExtension: \"Add a file extension\" TYPE [" + CURRENT + ": \"Keep current\", fa: \"FASTA\", gtf: \"GTF\"] DEFAULT " + CURRENT + " (The output file is named according to the last part of the URL. If it doesn't contain a correct file extension, select it here so that the file type is recognized correctly.)";
	}
	
	@Override
	protected void execute() { 
		updateStateToClient(JobState.RUNNING, "downloading");

		try {
			// file 
			File outputFile = new File(jobWorkDir, analysis.getOutputFiles().get(0).getFileName().getID()); 

			// parameters
			List<String> parameters = inputMessage.getParameters(JAVA_PARAMETER_SECURITY_POLICY, analysis);
			String urlString = parameters.get(0);
			String fileExtension = parameters.get(1);
			
			URL url = new URL(urlString);
			
			
			if (!allowedProtocols.contains(url.getProtocol())) {
				
				getResultMessage().setErrorMessage("CHIPSTER-NOTE: unsupported protocol: " + url.getProtocol());
				updateState(JobState.FAILED_USER_ERROR, "");
				return;
			}
			
			if (UrlTransferUtil.isLocalhost(url.getHost())) {

				getResultMessage().setErrorMessage("CHIPSTER-NOTE: not allowed to connect localhost: " + url.getHost());
				updateState(JobState.FAILED_USER_ERROR, "");
				return;
			}
			
			String datasetName = UrlTransferUtil.parseFilename(url);
			
			if (!CURRENT.equals(fileExtension)) {
				datasetName += "." + fileExtension;
			}

			FileUtils.copyURLToFile(url, outputFile, CONNECTION_TIMEOUT, READ_TIMEOUT);
			
			LinkedHashMap<String, String> nameMap = new LinkedHashMap<>();
			nameMap.put(outputFile.getName(), datasetName);
			
			JavaToolUtils.writeOutputDescription(jobWorkDir, nameMap);
						
		} catch (Exception e) {
			getResultMessage().setErrorMessage(Exceptions.getStackTrace(e));
			updateState(JobState.FAILED, "");
			return;
		}

		updateStateToClient(JobState.RUNNING, "download finished");
	}
}
