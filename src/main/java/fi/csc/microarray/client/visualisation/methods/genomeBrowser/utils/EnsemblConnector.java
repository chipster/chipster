package fi.csc.microarray.client.visualisation.methods.genomeBrowser.utils;

import java.sql.Connection;
import java.sql.DriverManager;

import com.mysql.jdbc.Statement;

public class EnsemblConnector {

	public static void main(String args[]){
		
		try {
			Statement stmt;

			//Register the JDBC driver for MySQL.
			Class.forName("com.mysql.jdbc.Driver");

			//Define URL of database server for
			// database named mysql on the localhost
			// with the default port number 3306.
			String url =
				"jdbc:mysql://ensembldb.ensembl.org:5306";

			//Get a connection to the database for a
			// user named root with a blank password.
			// This user is the default administrator
			// having full privileges to do anything.
			Connection con =
				DriverManager.getConnection(
						url,"anonymous", "");

			//Display URL and connection information
			System.out.println("URL: " + url);
			System.out.println("Connection: " + con);

			//Get a Statement object
			//stmt = con.createStatement();
			
			/*

# EnsEMBL KNOWN GENES - no coding sequence, just exons and introns, name(for example B3GALTL) and strand, this is most important

SELECT 
g.gene_id,
g.seq_region_start as gene_start,
g.seq_region_end as gene_end, 
g.seq_region_strand as strand, 
s.name as chr,
e.seq_region_start as exon_start, 
e.seq_region_end as exon_end
FROM gene g, seq_region s, exon_transcript et, exon e, coord_system c
WHERE 
g.status = 'KNOWN' AND
et.transcript_id = g.canonical_transcript_id AND
et.exon_id = e.exon_id AND
s.coord_system_id = c.coord_system_id AND
c.name='chromosome' AND
c.attrib = 'default_version'
ORDER BY s.name, gene_start, exon_start
LIMIT 0,1000
=cut


#GenScan PREDICTED GENES, left away for know

#EnsEMBL KNOWN TRANSCRIPTS - find and show UTR and CDS

SELECT exons.transcript_id, 
exons.transcript_start, 
exons.transcript_end, 
exons.strand, 
exons.exon_id, 
exons.exon_start, 
exons.exon_end, 
exons.chr,  
exons.symbol, 
tr.seq_start as CDS_start, 
tr.seq_end as CDS_end
 
FROM
(SELECT 
t.transcript_id, 
t.seq_region_start as transcript_start, 
t.seq_region_end as transcript_end, 
t.seq_region_strand as strand, 
e.exon_id,
e.seq_region_start as exon_start, 
e.seq_region_end as exon_end,
s.name as chr,
x.display_label as symbol

FROM transcript t, exon_transcript et, exon e, seq_region s, xref x, coord_system c
WHERE t.transcript_id = et.transcript_id AND
e.exon_id = et.exon_id AND
t.seq_region_id = s.seq_region_id AND
t.display_xref_id = x.xref_id AND
t.status = 'KNOWN' AND
c.name='chromosome' AND
c.attrib = 'default_version' AND
s.coord_system_id = c.coord_system_id AND
s.name<>'MT') as exons
LEFT JOIN translation tr
ON exons.transcript_id = tr.transcript_id AND 
exons.exon_id=tr.start_exon_id OR
exons.exon_id=tr.end_exon_id

ORDER BY transcript_start, exon_start, CDS_start
LIMIT 0,1000

# sno_miRNA

#CYTOBANDS

SELECT seq_region_start, seq_region_end, band, stain, name 
FROM homo_sapiens_core_56_37a.karyotype k
INNER JOIN homo_sapiens_core_56_37a.seq_region s
ON k.seq_region_id=s.seq_region_id;

			 */

			con.close();
		}catch( Exception e ) {
			e.printStackTrace();
		}//end catch
	}//end main
}//end class Jdbc11

