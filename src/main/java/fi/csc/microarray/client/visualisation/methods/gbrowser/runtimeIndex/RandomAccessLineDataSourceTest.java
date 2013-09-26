package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.DataUrl;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;

/**
 * Example output on ThinkPad X220:
 * 
 * true
First 10000 lines: 
	1: 1	processed_transcript	exon	11869	12227	.	+	.	 gene_id "ENSG00000223972"; transcript_id "ENST00000456328"; exon_number "1"; gene_name "DDX11L1"; gene_biotype "pseudogene"; transcript_name "DDX11L1-002"; exon_id "ENSE00002234944";
	10000: 1	protein_coding	CDS	6202188	6202388	.	-	0	 gene_id "ENSG00000116254"; transcript_id "ENST00000262450"; exon_number "15"; gene_name "CHD5"; gene_biotype "protein_coding"; transcript_name "CHD5-001"; protein_id "ENSP00000262450";
36 ms 
Seek to 100MB: 
	1: _name "WDR49"; gene_biotype "protein_coding"; transcript_name "WDR49-201"; exon_id "ENSE00002800706";
0 ms 
Seek to 300MB: 
	1: 00452240";
0 ms 
Last line: 
	1: Y	processed_pseudogene	exon	59001391	59001635	.	+	.	 gene_id "ENSG00000235857"; transcript_id "ENST00000431853"; exon_number "1"; gene_name "CTBP2P1"; gene_biotype "pseudogene"; transcript_name "CTBP2P1-001"; exon_id "ENSE00001794473";
0 ms 
First 10000 lines: 
	1: 1	processed_transcript	exon	11869	12227	.	+	.	 gene_id "ENSG00000223972"; transcript_id "ENST00000456328"; exon_number "1"; gene_name "DDX11L1"; gene_biotype "pseudogene"; transcript_name "DDX11L1-002"; exon_id "ENSE00002234944";
	10000: 1	nonsense_mediated_decay	exon	6453317	6453421	.	-	.	 gene_id "ENSG00000097021"; transcript_id "ENST00000418124"; exon_number "1"; gene_name "ACOT7"; gene_biotype "protein_coding"; transcript_name "ACOT7-003"; exon_id "ENSE00001556789";
2963 ms 
Seek to 100MB: 
	1: _name "WDR49"; gene_biotype "protein_coding"; transcript_name "WDR49-201"; exon_id "ENSE00002800706";
13 ms 
Seek to 300MB: 
	1: 00452240";
17 ms 
Last line: 
	1: Y	processed_pseudogene	exon	59001391	59001635	.	+	.	 gene_id "ENSG00000235857"; transcript_id "ENST00000431853"; exon_number "1"; gene_name "CTBP2P1"; gene_biotype "pseudogene"; transcript_name "CTBP2P1-001"; exon_id "ENSE00001794473";
6 ms 
Test done
 * 
 * @author klemela
 */
public class RandomAccessLineDataSourceTest {
	public static void main (String args[]) throws URISyntaxException, IOException, GBrowserException {
		
		//Functionality
		testFunctionality();
		
		
		//Performance
		URL fileUrl = new File(System.getProperty("user.home") + "/chipster/Homo_sapiens.GRCh37.69-sort.gtf").toURI().toURL();
		DataUrl fileDataUrl = new DataUrl(fileUrl, "file");
		
		RandomAccessLineDataSource file = new RandomAccessLineDataSource(fileDataUrl);		
		
		URL httpUrl = new URL("http://chipster-filebroker.csc.fi:7060/public/annotations/tmp/Homo_sapiens.GRCh37.69-sort.gtf");
		DataUrl httpDataUrl = new DataUrl(httpUrl, "http");
		RandomAccessLineDataSource http = new RandomAccessLineDataSource(httpDataUrl);
		
		manualTest(file);
		manualTest(http);
		
		iterativeTest(file, http);
	}

	private static void testFunctionality() throws URISyntaxException, IOException, GBrowserException {
				
		File testFile = RandomAccessLineReaderTest.getTestFile();
		DataUrl dataUrl = new DataUrl(testFile);
		List<String> lines = RandomAccessLineReaderTest.getTestReferenceList(testFile);
		
		RandomAccessLineDataSource dataSource = new RandomAccessLineDataSource(dataUrl);
		
		System.out.println(lines.get(lines.size() - 1).equals(dataSource.getLastLine()));
		
		testFile.delete();
	}

	private static void iterativeTest(RandomAccessLineDataSource file,
			RandomAccessLineDataSource http) throws IOException, GBrowserException {
		
		for (int i = 0; i < 1000; i++) {
			file.setLineReaderPosition(i);
			http.setLineReaderPosition(i);
			
			testCompare(file.getNextLine(), http.getNextLine(), "" + i);
		}
		
		for (int i = 1000000; i < 1001000; i++) {
			file.setLineReaderPosition(i);
			http.setLineReaderPosition(i);
			
			testCompare(file.getNextLine(), http.getNextLine(), "" + i);
		}
		
		testCompare(file.getLastLine(), http.getLastLine(), "lastLine");
		
		System.out.println("Test done");		
	}

	private static void testCompare(String line1, String line2, String description) {
		if (!line1.equals(line2)) {
			System.out.println("Lines are different!" + description);
			System.out.println("\t" + line1);
			System.out.println("\t" + line2);
		}
	}

	private static void manualTest(RandomAccessLineDataSource file)
			throws IOException, GBrowserException {
		long t = System.currentTimeMillis();
		
		System.out.println("First 10000 lines: ");
		file.setLineReaderPosition(0);
		System.out.println("\t1: " + file.getNextLine());
		
		for (int i = 0; i < 10000; i++) {
			file.getNextLine();
		}
		System.out.println("\t10000: " + file.getNextLine());
		
		System.out.println(System.currentTimeMillis() - t + " ms ");
		t = System.currentTimeMillis();
		
		System.out.println("Seek to 100MB: ");
		file.setLineReaderPosition(100*1024*1024);
		System.out.println("\t1: " + file.getNextLine());
		

		System.out.println(System.currentTimeMillis() - t + " ms ");
		t = System.currentTimeMillis();
		
		System.out.println("Seek to 300MB: ");
		file.setLineReaderPosition(300*1024*1024);
		System.out.println("\t1: " + file.getNextLine());

		System.out.println(System.currentTimeMillis() - t + " ms ");
		t = System.currentTimeMillis();
		
		System.out.println("Last line: ");
		System.out.println("\t1: " + file.getLastLine());
		
		System.out.println(System.currentTimeMillis() - t + " ms ");
		t = System.currentTimeMillis();
	}
}
