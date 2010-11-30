package fi.csc.microarray.gbrowser.index;

/**
 * Gene name search actions (creating database in memory, creating table,
 * updating it with gene information (chromosome, bp_start, bp_end, name),
 * getting gene coordinates by its name 
 * 
 * FIXME Check whether connections etc are evere released
 * 
 * @author zukauska
 */

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.ChunkDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AnnotationContents.Genome;

public class GeneIndexActions {

	private Connection conn;
	private Statement st;
	private static GeneIndexActions gia;
	private GetGeneIndexData getData;
	private static Genome selectedGenome;

	private GeneIndexActions(ChunkDataSource dataSource) throws ClassNotFoundException, SQLException {
		getData = new GetGeneIndexData(dataSource);
		Class.forName("org.h2.Driver");
		conn = DriverManager.getConnection("jdbc:h2:mem:GeneIndex", "sa", "");
		createTable();
		updateTable(getData.read());
	}

	/**
	 * creates database, table, updates table with gene indexes
	 * @throws SQLException 
	 * @throws ClassNotFoundException 
	 */
	public static GeneIndexActions getInstance(Genome genome, ChunkDataSource dataSource) throws ClassNotFoundException, SQLException {

		if (gia == null) {
			gia = new GeneIndexActions(dataSource);
			selectedGenome = genome;

		} else if (!selectedGenome.equals(genome)) {
			selectedGenome = genome;
			gia = new GeneIndexActions(dataSource);
			selectedGenome = genome;
		}

		return gia;
	}

	/**
	 * create table and index for geneName
	 * @throws SQLException 
	 */
	private void createTable() throws SQLException {

		st = conn.createStatement();
		st.execute("CREATE TABLE IF NOT EXISTS gene_name_index(" + "ID INT PRIMARY KEY auto_increment," + "chromosome VARCHAR(255),"
				+ "bp_start INT," + "bp_end INT," + "NAME VARCHAR(255),);"
				+ "create index IF NOT EXISTS geneName on gene_name_index(name);");
	}

	private void updateTable(List<RegionContent> indexList) throws SQLException {
		st = conn.createStatement();
		ResultSet rs = st.executeQuery("select * from gene_name_index limit 0,1");
		if (rs.next()) {
			st.execute("delete from gene_name_index");
		}
		st = conn.createStatement();
		for (RegionContent id : indexList) {
			st.executeUpdate("insert into gene_name_index values (null,'" + id.values.get(ColumnType.CHROMOSOME) + "',"
					+ id.values.get(ColumnType.BP_START) + "," + id.values.get(ColumnType.BP_END) + ",'"
					+ id.values.get(ColumnType.DESCRIPTION) + "' )");
		}
	}

	/**
	 * getting location of a gene
	 */
	public BpCoordRegion getLocation(String name, Chromosome chromosome) {
		try {
			st = conn.createStatement();
			ResultSet rs = st.executeQuery("select chromosome, bp_start,bp_end from gene_name_index " + "where name ='" + name
					+ "' and chromosome ='" + chromosome.toString() + "'");
			if (!rs.next()) {
				st = conn.createStatement();
				rs = st.executeQuery("select chromosome, bp_start,bp_end from gene_name_index " + "where name ='" + name + "'");
				rs.next();
			}
			return new BpCoordRegion(rs.getLong(2), rs.getLong(3), new Chromosome(rs.getString(1)));

		} catch (SQLException e) {
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

	public void closeConnection() {
		try {
			conn.close();
		} catch (SQLException e) {

		}
	}
}
