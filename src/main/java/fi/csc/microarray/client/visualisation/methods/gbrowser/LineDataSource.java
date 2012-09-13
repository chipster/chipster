package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.HttpURLConnection;
import java.net.URISyntaxException;
import java.net.URL;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;


/**
 * Simple data source for small non-binary files. These small files are stored in memory, so  only
 * simple initial line-by-line read through is needed.
 * 
 * @author Petri Klemelä
 *
 */
public class LineDataSource extends DataSource {
	
	public LineDataSource(URL url, Class<? extends AreaRequestHandler> requestHandler) throws FileNotFoundException, URISyntaxException {
		super(url, requestHandler);
	}

	BufferedReader reader;
	
	public String readLine() throws IOException {
		if (file != null) {

			if (reader == null) {
				reader = new BufferedReader(new FileReader(file));
			}


		} else {

			if (reader == null) {

				HttpURLConnection connection = (HttpURLConnection)url.openConnection();

				reader = new BufferedReader(new InputStreamReader(connection.getInputStream()));
			}
		}

		String result = reader.readLine();

		if (result == null) {
			reader.close();
			reader = null;
		}

		return result;
	}
}
