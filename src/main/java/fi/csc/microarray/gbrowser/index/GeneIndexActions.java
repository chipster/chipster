package fi.csc.microarray.gbrowser.index;

/**
 * Gene indexing tools. For a single genome, inserts gene names into 
 * an indexed database with their coordinates. Database is stored
 * in memory. Gene coordinates can be queried from the database via
 * this tool.
 * 
 * @author zukauska, Aleksi Kallio
 */
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.HashMap;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.ChunkDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AnnotationManager.Genome;

public class GeneIndexActions {

	private static HashMap<String, GeneIndexActions> instances = new HashMap<String, GeneIndexActions>();
	private PreparedStatement insertGeneStatement;
	private PreparedStatement selectStatement;

	protected GeneIndexActions(ChunkDataSource dataSource) throws ClassNotFoundException, SQLException {
		Class.forName("org.h2.Driver");
		Connection conn = DriverManager.getConnection("jdbc:h2:mem:GeneIndex", "", "");
		initialise(conn);
		
		if (dataSource != null) {
			GetGeneIndexData data = new GetGeneIndexData(dataSource);
			updateTable(data.read());
		}
	}

	private void initialise(Connection conn) throws SQLException {
		
		// Initialise the database 
		PreparedStatement createTableStatements = conn.prepareStatement(
				"DROP TABLE gene IF EXISTS;"
				+ "CREATE TABLE gene ("  
				+ "ID INT PRIMARY KEY AUTO_INCREMENT, chromosome VARCHAR(255),"
				+ "bp_start INT, bp_end INT, name VARCHAR(255));"
				+ "DROP INDEX gene_name_index IF EXISTS;"
				+ "CREATE INDEX gene_name_index ON gene(name);");

		createTableStatements.execute();

		// Initialise statement for inserting genes
		this.insertGeneStatement = conn.prepareStatement("INSERT INTO gene VALUES (NULL, ?, ?, ?, ?);");

		// Initialise statement for finding gene locations
		this.selectStatement = conn.prepareStatement("SELECT chromosome, bp_start, bp_end FROM gene WHERE name = ?");

	}

	/**
	 * Returns gene index for given genome. If index has not yet been initialised, initialises it with
	 * given dataSource.
	 */
	public static GeneIndexActions getInstance(Genome genome, ChunkDataSource dataSource) throws ClassNotFoundException, SQLException {

		String genomeString = genome.toString();
		
		if (!instances.containsKey(genomeString)) {
			instances.put(genomeString, new GeneIndexActions(dataSource));
		}

		return instances.get(genomeString);
	}

	private void updateTable(List<RegionContent> indexList) throws SQLException {
		for (RegionContent id : indexList) {
			String chromosome = id.values.get(ColumnType.CHROMOSOME).toString();
			String bpStart = id.values.get(ColumnType.BP_START).toString();
			String bpEnd = id.values.get(ColumnType.BP_END).toString();
			String name = id.values.get(ColumnType.DESCRIPTION).toString();
			insertGene(chromosome, bpStart, bpEnd, name);
		}
	}

	protected void insertGene(String chromosome, String bpStart, String bpEnd, String name) throws SQLException {
		insertGeneStatement.setString(1, chromosome);
		insertGeneStatement.setString(2, bpStart);
		insertGeneStatement.setString(3, bpEnd);
		insertGeneStatement.setString(4, name.toUpperCase());
		insertGeneStatement.executeUpdate();
	}

	/**
	 * getting location of a gene
	 * @throws SQLException 
	 */
	public BpCoordRegion getLocation(String name) throws SQLException {
		selectStatement.setString(1, name.toUpperCase());
		ResultSet rs = selectStatement.executeQuery();
		if (rs.next()) {
			return new BpCoordRegion(rs.getLong(2), rs.getLong(3), new Chromosome(rs.getString(1)));
		} else {
			return null;
		}
	}

	public static boolean checkIfNumber(String name) {
		try {
			Integer.parseInt(name);
			return true;
		} catch (NumberFormatException e) {
			try {
				Long.parseLong(name);
				return true;
			} catch (NumberFormatException e1) {
				return false;
			}
		}
	}
}
