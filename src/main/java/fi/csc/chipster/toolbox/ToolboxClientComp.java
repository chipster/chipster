package fi.csc.chipster.toolbox;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.net.URL;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

import javax.ws.rs.NotFoundException;
import javax.ws.rs.client.Client;
import javax.ws.rs.client.ClientBuilder;
import javax.ws.rs.client.WebTarget;
import javax.ws.rs.core.MediaType;

import org.apache.log4j.Logger;

import com.fasterxml.jackson.core.JsonParseException;
import com.fasterxml.jackson.databind.JsonMappingException;

public class ToolboxClientComp {

	private String baseUri;
	private Client client;

	private final static String MODULES_ZIP_PATH = "/modules/zip";
	
	private static final Logger logger = Logger.getLogger(ToolboxClientComp.class);
	
	
	public ToolboxClientComp(String toolboxUri) {
		this.baseUri = toolboxUri;
		this.client = ClientBuilder.newClient();
	}

	public ToolboxTool getTool(String toolId) throws IOException {

		WebTarget serviceTarget = client.target(baseUri).path("tools/" + toolId);

		String json; 
		try {
			json = serviceTarget.request(MediaType.APPLICATION_JSON).get(String.class);
		} catch (NotFoundException nfe) {
			return null;
		}

		ToolboxTool tool = ToolboxRestUtils.parseJson(ToolboxTool.class, json, false);

		return tool;
	}

	public void close() {
		client.close();
	}

	
	public void getToolboxModules(File jobToolboxDir) throws IOException {
		unzip(baseUri + MODULES_ZIP_PATH, jobToolboxDir);
		// FIXME fix permissions
	}
	
	
	private void unzip(String zipFilePath, File destDirectory) throws IOException {
		long startTime = System.currentTimeMillis();

		File destDir = destDirectory;
		if (!destDir.exists()) {
			destDir.mkdir();
		}

		URL url = new URL(zipFilePath);

		try (ZipInputStream zipIn = new ZipInputStream(new BufferedInputStream(url.openStream(), 1024))) {
			ZipEntry entry = zipIn.getNextEntry();

			// iterates over entries in the zip file
			while (entry != null) {
				String filePath = destDirectory + File.separator + entry.getName();
				if (!entry.isDirectory()) {
					// if the entry is a file, extracts it
					extractFile(zipIn, filePath);
				} else {
					// if the entry is a directory, make the directory
					File dir = new File(filePath);
					dir.mkdir();
				}
				zipIn.closeEntry();
				entry = zipIn.getNextEntry();
			}
			zipIn.close();
		}
		
		logger.info("toolbox download took " + (System.currentTimeMillis() - startTime) + " ms");
	}

	/**
	 * Extracts a zip entry (file entry)
	 * @param zipIn
	 * @param filePath
	 * @throws IOException
	 */
	private void extractFile(ZipInputStream zipIn, String filePath) throws IOException {
		try (BufferedOutputStream bos = new BufferedOutputStream(new FileOutputStream(filePath))) {
			byte[] bytesIn = new byte[4096];
			int read = 0;
			while ((read = zipIn.read(bytesIn)) != -1) {
				bos.write(bytesIn, 0, read);
			}
		}
	}
	
	
	
	public static void main(String args[]) throws JsonParseException, JsonMappingException, IOException {
		
		ToolboxClientComp toolboxClient = new ToolboxClientComp("http://localhost:8008/toolbox");
		try {
			ToolboxTool tool = toolboxClient.getTool("norm-affy.R");
		
			System.out.println(tool.getSadlString());
		} finally {
			toolboxClient.close();
		}
	}

}
