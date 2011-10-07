package fi.csc.microarray.client.session;

import java.io.File;
import java.io.IOException;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

import javax.xml.XMLConstants;
import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.transform.stream.StreamSource;
import javax.xml.validation.Schema;
import javax.xml.validation.SchemaFactory;
import javax.xml.validation.Validator;

import org.apache.log4j.Logger;
import org.xml.sax.SAXException;


public class UserSession {

	public static final String SESSION_FILE_EXTENSION = "zip";
	public static final String SESSION_DATA_FILENAME = "session.xml";
	
	protected static final String SESSION_BACKUP_PREFIX = "backup_session";
	
	private static final Logger logger = Logger.getLogger(UserSession.class);
	public static final String ROOT_FOLDER_ID = "0";
	public static final Integer SESSION_VERSION = 1;
	
	
	public static boolean isValidSessionFile(File file) {
		if (file == null) {
			return false;
		}
		
		// does the file exist?
		if (!file.exists()) {
			return false;
		}

		// correct extension?
		if (!file.getName().endsWith("." + SESSION_FILE_EXTENSION)) {
			return false;
		}

		// is it a zip?
		ZipFile zipFile = null;
		try {
			try {
				 zipFile= new ZipFile(file);
			} catch (ZipException e) {
				return false;
			} catch (IOException e) {
				return false;
			}

			// does it contain the session metadata file
			if (zipFile.getEntry(SESSION_DATA_FILENAME) == null) {
				return false;
			}
		
		} finally {
			if (zipFile != null) {
				try {
					zipFile.close();
				} catch (IOException e) {
					logger.warn("could not close zip file: " + file.getName());
				}
			}
		}

		return true;
	}
	
	
    public static boolean validateMetadataFile() throws IOException, SAXException  {

        // 1. Specify you want a factory for RELAX NG
        SchemaFactory factory = SchemaFactory.newInstance(XMLConstants.W3C_XML_SCHEMA_NS_URI);
        
        // 2. Load the specific schema you want. 
        // Here I load it from a java.io.File, but we could also use a 
        // java.net.URL or a javax.xml.transform.Source
//        File schemaLocation = new File("/opt/xml/docbook/rng/docbook.rng");
        
        // 3. Compile the schema.
        Schema schema = factory.newSchema(new StreamSource(UserSession.class.getResourceAsStream("/session.xsd")));
    
        // 4. Get a validator from the schema.
        Validator validator = schema.newValidator();
        
        // 5. Parse the document you want to check.
        
        // 6. Check the document
        try {
        	validator.validate(new StreamSource(UserSession.class.getResourceAsStream("/session.xml")));
            System.out.println("input is valid.");
        }
        catch (SAXException ex) {
            System.out.println("input is not valid because ");
            System.out.println(ex.getMessage());
            return false;
        }  
        return true;
        
    }
	
    
    public static void main(String[] args) throws IOException, SAXException {
		validateMetadataFile();
	}
	public static JAXBContext getJAXBContext() throws JAXBException {
		return JAXBContext.newInstance("fi.csc.microarray.client.session.schema");
	}


	public static Schema getSchema() throws SAXException {
        SchemaFactory factory = SchemaFactory.newInstance(XMLConstants.W3C_XML_SCHEMA_NS_URI);
		return factory.newSchema(new StreamSource(UserSession.class.getResourceAsStream("session.xsd")));
	}


	public static File findBackupFile(File directory, boolean oneAfterLast) {
		
		File sessionFile = null;
		File candidateSessionFile = new File(directory, SESSION_BACKUP_PREFIX + "." + SESSION_FILE_EXTENSION);

		for (int i = 0; candidateSessionFile.exists(); i++) {
			sessionFile = candidateSessionFile;
			candidateSessionFile = new File(directory, SESSION_BACKUP_PREFIX + i + "." + SESSION_FILE_EXTENSION);
		}
		
		if (oneAfterLast) {
			sessionFile = candidateSessionFile;
		}
		return sessionFile;
	}
}
