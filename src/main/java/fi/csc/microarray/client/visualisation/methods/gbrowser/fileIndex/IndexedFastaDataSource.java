package fi.csc.microarray.client.visualisation.methods.gbrowser.fileIndex;

import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;

import net.sf.picard.PicardException;
import net.sf.picard.reference.ChipsterIndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.ByteDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.LineDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.ChromosomeNameUnnormaliser;

/**
 * 
 * @author Petri Klemelä
 *
 */
public class IndexedFastaDataSource extends DataSource {
	
	private ChromosomeNameUnnormaliser chromosomeNameUnnormaliser = ChromosomeNameUnnormaliser.newIdentityPreversingUnnormaliser();

	private URL dataUrl;
    private URL indexUrl;
	private ChipsterIndexedFastaSequenceFile picard;

	public IndexedFastaDataSource(URL data, URL index) throws URISyntaxException, IOException {
        super(data, null);
        
        this.dataUrl = data;
        this.indexUrl = index;
        
		ByteDataSource dataSource = new ByteDataSource(dataUrl);
		LineDataSource indexSource = new LineDataSource(indexUrl);

		picard = new ChipsterIndexedFastaSequenceFile(dataSource, indexSource);
		
		// Create unnormaliser for this naming convention
		this.chromosomeNameUnnormaliser = new ChromosomeNameUnnormaliser(picard.getContigs());		
    }
	
	public URL getIndex() {
		return indexUrl;
	}

	public URL getUrl() {
		return dataUrl;
	}

	public String query(Chromosome chr, Long start, Long end) {
		
		String chrString = chromosomeNameUnnormaliser.unnormalise(chr);
		
		try {
			ReferenceSequence picardSequence = picard.getSubsequenceAt(chrString, start, end);
			return new String(picardSequence.getBases());

		} catch (PicardException e) {				
			e.printStackTrace(); //Catch "Query asks for data past end of contig" to prevent this thread from ending
			return null;
		}		
	}
}
