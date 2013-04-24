package fi.csc.microarray.client.visualisation.methods.gbrowser.stack;

import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;

import org.broad.tribble.readers.TabixReader;

import net.sf.picard.PicardException;
import net.sf.picard.reference.ReferenceSequence;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.BpCoordFileRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.IndexedFastaDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.TabixDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Exon;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Gene;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.GeneRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.GeneSet;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.ParsedFileResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.TabixUtil;

public class GtfTabixToFeatureConversion extends SingleThreadAreaRequestHandler {

	private TabixDataSource dataSource;
	private GtfLineParser parser = new GtfLineParser();

	public GtfTabixToFeatureConversion(URL data, URL index, final GBrowser browser) {

		super(null, null);

		try {			
			this.dataSource = new TabixDataSource(data, index);

		} catch (URISyntaxException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	@Override
	public void clean() {
	}


	@Override
	protected void processAreaRequest(AreaRequest request) {

		super.processAreaRequest(request);

		if (request.getStatus().poison) {
			return;
		}
		
		List<RegionContent> resultList = null;

		if (request instanceof GeneRequest) {

			try {
				resultList = processGeneSearch((GeneRequest)request);
				
			} catch (IOException e) {
				e.printStackTrace();
			}

		} else { 

			resultList = new LinkedList<RegionContent>();

			GeneSet genes;
			try {
				genes = fetchExons(request);

				for (Gene gene : genes.values()) {

					LinkedHashMap<ColumnType, Object> values = new LinkedHashMap<ColumnType, Object>();

					values.put(ColumnType.VALUE, gene);

					resultList.add(new RegionContent(gene.getRegion(), values));
				}
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		createAreaResult(new AreaResult(request.getStatus(), resultList));	
	}
	
	private List<RegionContent> processGeneSearch(GeneRequest areaRequest) throws IOException {

		String searchString = areaRequest.getSearchString().toLowerCase();
		Chromosome chr = areaRequest.start.chr;


		Region region = new Region(1l, Long.MAX_VALUE, chr);

		GeneSet genes = fetchExons(region);
		
		List<RegionContent> resultList = new LinkedList<RegionContent>();

		for (Gene gene : genes.values()) {

			if (gene.getName() != null && gene.getName().toLowerCase().equals(searchString)) {

				LinkedHashMap<ColumnType, Object> values = new LinkedHashMap<ColumnType, Object>();

				values.put(ColumnType.VALUE, gene);
				resultList.add(new RegionContent(gene.getRegion(), values));
			}
		}

		return resultList;
	}

	/**
	 * Find exons in a given range.
	 * 
	 * @param request
	 * @return
	 * @throws IOException 
	 */
	public GeneSet fetchExons(Region request) throws IOException {
		// Read the given region
		TabixReader.Iterator iter = TabixUtil.getTabixIterator(dataSource, request);

		GeneSet genes = new GeneSet();		
		String line;
		
		if (iter != null) { //null if there isn't such chromosome in annotations

			while ((line = iter.next()) != null) {
				
				parser.setLine(line);
				
				
				parseGtfLine(line, genes);
			}
		}

		return genes;
	}

	private void parseGtfLine(String line, GeneSet genes) {
		
		String[] cols;

		String[] ids;

		cols = line.split("\t");

		String chr = cols[0];
		String biotype = cols[1];
		String feature = cols[2];
		String exonStart = cols[3];
		String exonEnd = cols[4];
		String strand = cols[6];

		ids = parseIds(cols[8]);

		String geneId = ids[0];
		String transcId = ids[1];
		String exonIndex = ids[2];
		String geneName = ids[3];
		String transcName = ids[4];
		
		if ("exon".equals(feature) || "CDS".equals(feature)) {

			Region region = new Region(Long.parseLong(exonStart), Long.parseLong(exonEnd), 
					new Chromosome(chr), getStrand(strand));

			Exon exon = new Exon(region, feature, Integer.parseInt(exonIndex));

			genes.addExon(exon, geneId, transcId, geneName, transcName, biotype);
			
		} else 	if (feature.startsWith("GenBank")) {

			Region region = new Region(Long.parseLong(exonStart), Long.parseLong(exonEnd), 
					new Chromosome(chr), getStrand(strand));
			
			if (geneId == null || transcId == null) {
				return;
			}
			
			if ("GenBank gene".equals(feature)) {
				feature = "exon";
			} else if ("GenBank CDS".equals(feature)) {
				feature = "CDS";
			} else {
				geneId = feature + geneId;
				transcId = feature + transcId;
				
				if (geneName != null) {
					geneName = feature + " " + geneName;
				}
				
				if (transcName != null) {
					transcName = feature + " " + transcName;
				}
				
				feature = "exon";
			}
			
			exonIndex = "1";
			

			Exon exon = new Exon(region, feature, Integer.parseInt(exonIndex));

			genes.addExon(exon, geneId, transcId, geneName, transcName, biotype);
		}
	} 

	private static Strand getStrand(String strand) {

		if ("+".equals(strand)) {
			return Strand.FORWARD;
		} else if ("-".equals(strand)) {
			return Strand.REVERSE;
		} else if (".".equals(strand)) {
			return Strand.NONE;
		}
		return Strand.UNRECOGNIZED;
	}

	private static final String[] ID_FIELDS = { "gene_id", "transcript_id", "exon_number", "gene_name", "transcript_name" };

	public static String[] parseIds(String ids) {

		String[] split = ids.split(";");
		String[] result = new String[ID_FIELDS.length];

		String key = null;
		String value = null;
		int indexOfQuotationMark = 0;

		for (int i = 0; i < split.length; i++) {

			indexOfQuotationMark = split[i].indexOf("\"");

			key = split[i].substring(1, indexOfQuotationMark - 1);
			value = split[i].substring(indexOfQuotationMark + 1, split[i]
					.lastIndexOf("\""));

			for (int fieldNumber = 0; fieldNumber < ID_FIELDS.length; fieldNumber++) {
				if (ID_FIELDS[fieldNumber].equals(key)) {
					result[fieldNumber] = value;
				}
			}
		}

		return result;
	}

	public class Transcript {

		Region region;
		List<Exon> exons;
	}

	public String toString() {
		return this.getClass().getName() + " - " + dataSource;
	}
}
