package fi.csc.microarray.client.visualisation.methods.gbrowser.stack;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URISyntaxException;
import java.net.URL;

import fi.csc.microarray.util.IOUtils;

public class FileLineReader implements LineReader {

	private BufferedReader lineReader;
	private File file;

	public FileLineReader(URL url) throws FileNotFoundException, URISyntaxException {
		this.file = new File(url.toURI());
	}	
	
	public void setPosition(long position) throws IOException {
		
		if (lineReader != null) {
			IOUtils.closeIfPossible(lineReader);
		}

		FileInputStream in = new FileInputStream(file);
			in.skip(position);
			lineReader = new BufferedReader(new InputStreamReader(in));
	}


	public void close() {
		if (lineReader != null) {
			try {
				lineReader.close();
			} catch (IOException e) {
				//No problem
			}
			lineReader = null;
		}
	}
	
	public long length() throws IOException {
		return file.length();
	}

	@Override
	public String readLine() throws IOException {
		return lineReader.readLine();
	}
}
