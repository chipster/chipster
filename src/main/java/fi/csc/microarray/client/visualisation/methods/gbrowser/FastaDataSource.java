package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;

/**
 * Data source for indexed SAM compatible files (.sam/.bam). Because index is required, in practice
 * supports only .bam.
 * 
 * @author Petri Klemelä
 *
 */
public class FastaDataSource extends DataSource {
	

	private Map<Chromosome, ChunkDataSource> sources = new TreeMap<Chromosome, ChunkDataSource>();
	

    public FastaDataSource() throws FileNotFoundException {
        super(new File(""));
    }

	public Set<Entry<Chromosome, ChunkDataSource>> entrySet() {
		return sources.entrySet();
	}

	public void put(Chromosome chromosome, File file) {
		try {
			sources.put(chromosome, new ChunkDataSource(file, null));
		} catch (FileNotFoundException e) {
			
			e.printStackTrace();
		}
	}
}
