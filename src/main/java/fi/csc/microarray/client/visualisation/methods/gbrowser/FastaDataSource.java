package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.io.FileNotFoundException;
import java.net.MalformedURLException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.FastaHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;

/**
 * Data source for fasta files. Each chromosome has a separate file that are stored in Map 
 * data structure. ChunkDataSources are used to get the random access functionality.
 * Having multiple DataSources inside one DataSource doesn't fit especially nicely with the
 * original concept of DataSource, but was needed to disguise this group of files into a form
 * that is allowed to travel through other layers of the system.
 * 
 * @author Petri Klemelä
 *
 */
public class FastaDataSource extends DataSource {
	

	private Map<Chromosome, ChunkDataSource> sources = new TreeMap<Chromosome, ChunkDataSource>();
	

    public FastaDataSource() throws FileNotFoundException, MalformedURLException, URISyntaxException {
        super(null, FastaHandlerThread.class);
    }

	public Set<Entry<Chromosome, ChunkDataSource>> entrySet() {
		return sources.entrySet();
	}

	public void put(Chromosome chromosome, URL url) throws URISyntaxException {
		try {
			sources.put(chromosome, new ChunkDataSource(url, null, FastaHandlerThread.class));
		} catch (FileNotFoundException e) {
			
			e.printStackTrace();
		}
	}
}
